from __future__ import absolute_import
from __future__ import print_function

import multiprocessing
import sys, time
from xfab import tools, detector
from fabio import edfimage, tifimage
import gzip  # to write .edf.gz
from . import variables, check_input
from scipy import ndimage
from scipy.stats import norm
import numpy as n
from PIL import Image

A_id = variables.refarray().A_id

# Set for 2kx2k detector, 10 deg corresponding to 1k
# and then reduced by a factor of 2-3 to get smooth peaks
Delta_tth = 0.005  # approx size of one pixel in deg tth
Delta_eta = 0.02  # approx size of one pixel in deg eta


class make_image:
    def __init__(self, graindata, killfile=None, frame_number=None):
        self.graindata = graindata
        self.killfile = killfile
        peakshape = self.graindata.param['peakshape']
        if peakshape[0] == 0:  # spike peak, 2x2 pixels
            self.peak_add = 1
            self.frame_add = 1
            self.peakwsig = 0
        elif peakshape[0] == 1:  # 3d Gaussian peak
            self.peak_add = max(1, int(round(peakshape[1])))
            self.frame_add = max(1, int(round(peakshape[1])))
            self.peakwsig = peakshape[2]
        elif peakshape[0] == 3:  # 3d Gaussian peak in 2theta,eta,omega
            self.peak_add = 1
            self.frame_add = 1
            self.cen_tth = int(1.5 * peakshape[1] / Delta_tth)
            self.frame_tth = 2 * self.cen_tth + 1
            self.fwhm_tth = peakshape[1] / Delta_tth
            self.cen_eta = int(1.5 * peakshape[2] / Delta_eta)
            self.frame_eta = 2 * self.cen_eta + 1
            self.fwhm_eta = peakshape[2] / Delta_eta
            self.aw_tth_eta = n.zeros((self.frame_tth, self.frame_eta))
            self.raw_tth_eta[self.cen_tth, self.cen_eta] = 1
            self.filter_tth_eta = ndimage.gaussian_filter(self.raw_tth_eta, [0.5 * self.fwhm_tth, 0.5 * self.fwhm_eta])
            self.peakwsig = 1.

        self.framedimy = int(self.graindata.param['dety_size'] + 2 * self.frame_add)
        self.framedimz = int(self.graindata.param['detz_size'] + 2 * self.frame_add)

        self.totalrefl = 0
        if frame_number is None:
            self.no_frames = list(range(len(self.graindata.frameinfo)))
            print('Generating diffraction images')
        else:
            self.no_frames = [frame_number]

    def make_one_frame(self, this_frame, is_multi):
        peakshape = self.graindata.param['peakshape']
        check_input.interrupt(self.killfile)
        t1 = time.clock()
        nrefl = 0
        frame = n.zeros((self.framedimy, self.framedimz))
        omega = self.graindata.frameinfo[this_frame].omega
        omega_step = self.graindata.param['omega_step']
        # Jettes hack to add relative movement of sample and detector, modelled to be Gaussian in y and z direction with a spread of 1 micron
        # movement of 1 micron along x judged to be irrelevant, at least for farfield data
        y_move = n.random.normal(0, 1. / self.graindata.param['dety_size'])
        z_move = n.random.normal(0, 1. / self.graindata.param['detz_size'])
        # loop over grains
        for j in range(self.graindata.param['no_grains']):
            # loop over reflections for each grain
            gr_pos = n.array(self.graindata.param['pos_grains_%s' % j])
            for k in range(len(self.graindata.grain[j].refs)):
                # exploit that the reflection list is sorted according to omega
                if self.graindata.grain[j].refs[k, A_id['omega']] * 180 / n.pi > \
                        omega + omega_step + 2 * self.peakwsig:
                    break
                elif self.graindata.grain[j].refs[k, A_id['omega']] * 180 / n.pi < \
                        omega - 2 * self.peakwsig:
                    continue
                dety = self.graindata.grain[j].refs[k, A_id['detyd']]  # must be spot position after
                detz = self.graindata.grain[j].refs[k, A_id['detzd']]  # applying spatial distortion
                # apply hack
                #                   dety = self.graindata.grain[j].refs[k,A_id['dety']] + y_move
                #                   detz = self.graindata.grain[j].refs[k,A_id['detz']] + z_move
                ndety = int(round(dety))
                ndetz = int(round(detz))
                yrange = list(range(ndety + self.frame_add - self.peak_add, ndety + self.frame_add + self.peak_add + 1))
                zrange = list(range(ndetz + self.frame_add - self.peak_add, ndetz + self.frame_add + self.peak_add + 1))
                intensity = int(round(self.graindata.grain[j].refs[k, A_id['Int']]))
                nrefl = nrefl + 1
                self.totalrefl = self.totalrefl + 1
                # Gaussian along omega
                if peakshape[0] == 1 or peakshape[0] == 3:
                    fraction = norm.cdf(
                        (omega - self.graindata.grain[j].refs[k, A_id['omega']] * 180 / n.pi + omega_step) / (
                                0.5 * self.peakwsig)) \
                               - norm.cdf(
                        (omega - self.graindata.grain[j].refs[k, A_id['omega']] * 180 / n.pi) / (0.5 * self.peakwsig))
                else:
                    fraction = 1.
                if peakshape[0] == 3:
                    # Gaussian peaks along 2theta,eta
                    tth = self.graindata.grain[j].refs[k, A_id['tth']]
                    eta = self.graindata.grain[j].refs[k, A_id['eta']]
                    Om = tools.form_omega_mat_general(self.graindata.grain[j].refs[k, A_id['omega']], 0,
                                                      -1. * self.graindata.param['wedge'] * n.pi / 180.)
                    [tx, ty, tz] = n.dot(Om, gr_pos)
                    for t in range(self.frame_tth):
                        tth_present = tth + (t - self.cen_tth) * Delta_tth * n.pi / 180.
                        for e in range(self.frame_eta):
                            eta_present = eta + (e - self.cen_eta) * Delta_eta * n.pi / 180.
                            [dety_present, detz_present] = detector.det_coor2(tth_present,
                                                                              eta_present,
                                                                              self.graindata.param['distance'],
                                                                              self.graindata.param['y_size'],
                                                                              self.graindata.param['z_size'],
                                                                              self.graindata.param['dety_center'],
                                                                              self.graindata.param['detz_center'],
                                                                              self.graindata.R,
                                                                              tx,
                                                                              ty,
                                                                              tz, )

                            if self.graindata.param['spatial'] is not None:
                                from ImageD11 import blobcorrector
                                self.spatial = blobcorrector.correctorclass(self.graindata.param['spatial'])
                                # To match the coordinate system of the spline file
                                # SPLINE(i,j): i = detz; j = (dety_size-1)-dety
                                # Well at least if the spline file is for frelon2k
                                (x, y) = detector.detyz_to_xy([dety_present, detz_present],
                                                              self.graindata.param['o11'],
                                                              self.graindata.param['o12'],
                                                              self.graindata.param['o21'],
                                                              self.graindata.param['o22'],
                                                              self.graindata.param['dety_size'],
                                                              self.graindata.param['detz_size'])
                                # Do the spatial distortion
                                (xd, yd) = self.spatial.distort(x, y)

                                # transform coordinates back to dety,detz
                                (dety_present, detz_present) = detector.xy_to_detyz([xd, yd],
                                                                                    self.graindata.param['o11'],
                                                                                    self.graindata.param['o12'],
                                                                                    self.graindata.param['o21'],
                                                                                    self.graindata.param['o22'],
                                                                                    self.graindata.param['dety_size'],
                                                                                    self.graindata.param['detz_size'])

                            y = int(round(dety_present))
                            z = int(round(detz_present))
                            try:
                                frame[y + self.frame_add, z + self.frame_add] = frame[
                                                                                    y + self.frame_add, z + self.frame_add] + fraction * intensity * \
                                                                                self.filter_tth_eta[t, e]
                            except:
                                # FIXME                                    print("Unhandled exception in make_image.py")
                                pass
                else:
                    # Generate spikes, 2x2 pixels
                    for y in yrange:
                        for z in zrange:
                            if 0 < y < self.framedimy and 0 < z < self.framedimz and abs(
                                    dety + self.frame_add - y) < 1 and abs(detz + self.frame_add - z) < 1:
                                #                                   frame[y-1,z] = frame[y-1,z] + fraction*intensity*(1-abs(dety+frame_add-y))*(1-abs(detz+frame_add-z))
                                y = int(round(y))
                                z = int(round(z))
                                frame[y, z] = frame[y, z] + fraction * intensity * (
                                        1 - abs(dety + self.frame_add - y)) * (
                                                      1 - abs(detz + self.frame_add - z))

        # 2D Gaussian on detector
        if peakshape[0] == 1:
            frame = ndimage.gaussian_filter(frame, peakshape[1] * 0.5)
        # add background
        if self.graindata.param['bg'] > 0:
            frame = frame + self.graindata.param['bg'] * n.ones((self.framedimy, self.framedimz))
        # add noise
        if self.graindata.param['noise'] != 0:
            frame = n.random.poisson(frame)
        # apply psf
        if self.graindata.param['psf'] != 0:
            frame = ndimage.gaussian_filter(frame, self.graindata.param['psf'] * 0.5)
        # resize, convert to integers and flip to same orientation as experimental frames
        frame = frame[self.frame_add:self.framedimy - self.frame_add, self.frame_add:self.framedimz - self.frame_add]

        # limit values above 16 bit to be 16bit
        frame = n.clip(frame, 0, 2 ** 16 - 1)
        # convert to integers
        frame = n.uint16(frame)

        # flip detector orientation according to input: o11, o12, o21, o22
        frame = detector.trans_orientation(frame,
                                           self.graindata.param['o11'],
                                           self.graindata.param['o12'],
                                           self.graindata.param['o21'],
                                           self.graindata.param['o22'],
                                           'inverse')
        # Output frames
        if '.edf' in self.graindata.param['output']:
            self.write_edf(this_frame, frame)
        if '.edf.gz' in self.graindata.param['output']:
            self.write_edf(this_frame, frame, usegzip=True)
        if '.tif' in self.graindata.param['output']:
            self.write_tif(this_frame, frame)
        if '.tif16bit' in self.graindata.param['output']:
            self.write_tif16bit(this_frame, frame)
        if not is_multi:
            print('\rDone frame %i took %8f s' % (this_frame + 1, time.clock() - t1), end=' ')
        sys.stdout.flush()

    def make_images(self):
        """
        makeimage script produces edf diffraction images using the reflection information
        
                Henning Osholm Sorensen, June 23, 2006.
        python translation Jette Oddershede, March 31, 2008
        """
        make_image_start = time.time()
        for i in self.no_frames:
            self.make_one_frame(i, is_multi=False)
        make_image_end = time.time()
        time_taken = make_image_end - make_image_start
        print("\nFinished making images. Took", time_taken, "seconds")

    def make_images_multi(self):
        """
        Multithreaded version of self.make_images
        """
        from multiprocessing import Pool
        from functools import partial
        from contextlib import contextmanager

        @contextmanager
        def poolcontext(*args, **kwargs):
            pool = multiprocessing.Pool(*args, **kwargs)
            yield pool
            pool.terminate()

        make_image_start = time.time()
        with poolcontext(processes=multiprocessing.cpu_count() - 1) as pool:
            results = pool.map(partial(self.make_one_frame, is_multi=True), self.no_frames)
        make_image_end = time.time()
        time_taken = make_image_end - make_image_start
        print("Finished making images. Took", time_taken, "seconds")

    def write_edf(self, framenumber, frame, usegzip=False):
        e = edfimage.edfimage()
        e.data = frame
        edim2, edim1 = frame.shape
        e.header = {}
        e.header['origin'] = 'PolyXSim'
        e.header['Dim_1'] = edim1
        e.header['Dim_2'] = edim2
        e.header['col_end'] = edim1 - 1
        e.header['row_end'] = edim2 - 1
        e.header['DataType'] = 'UnsignedShort'
        e.header['Image'] = 1
        e.header['ByteOrder'] = 'Low'
        e.header['time'] = time.asctime()
        e.header['Omega'] = self.graindata.frameinfo[framenumber].omega + \
                            self.graindata.param['omega_step'] / 2.0
        e.header['OmegaStep'] = self.graindata.param['omega_step']
        e.header['grainfile'] = '%s/%s_%0.4dgrains.txt' \
                                % (self.graindata.param['direc'], self.graindata.param['stem'],
                                   self.graindata.param['no_grains'])
        fname = '%s%s' % (self.graindata.frameinfo[framenumber].name, '.edf')
        if usegzip:
            fobj = gzip.GzipFile(fname + ".gz", "wb")
            e.write(fobj)
            fobj.close()
        else:
            e.write(fname)

    def write_tif(self, framenumber, frame):
        e = tifimage.tifimage()
        e.data = frame
        e.write('%s%s' % (self.graindata.frameinfo[framenumber].name, '.tif'))

    def write_tif16bit(self, framenumber, frame):
        size = frame.shape[:2][::-1]
        pilimage = Image.frombuffer('I', size, frame.tostring(), "raw", 'I;16', 0, 1)
        pilimage.save('%s%s' % (self.graindata.frameinfo[framenumber].name, '.tif'))
