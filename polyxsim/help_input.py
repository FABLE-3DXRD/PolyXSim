def show_input():
    """
    This function returns information about the structure, keywords and syntax
    of the PolyXSim input file as one string.
    """
    return \
    """
The input to PolyXSim is given as an ascii file. The file can be given any name.

A detailed definition of the experimental setup and parameters can be found on the wiki http://fable.sourceforge.net.

1.0 INPUT SYNTAX

The input is given by a set of input keywords followed by one or several values dependent on the type of keyword specified.
The ordering of the keywords is irrelevant.

The syntax is as follows
  * keyword ''values'' [units]

Followed by an input example:

keyword value                    # some explanation
if the value is a string in should be given in quotation marks 'value'

keyword 'value'                  # string values should en enclosed by quotation marks

if the value is a numeric number the format is free. [[br]]
E.g.
keyword 10
keyword 10.0
keyword 1e1
keyword 10**1
keyword 20/2.0

anything goes - or almost anything!

2.0   INPUT PARAMETERS

2.1   Instrumentation 

2.1.1 Beam specs 
  * wavelength ''wavelength'' [AA]
  * beam_width ''width-of-beam'' [mm]
  * beamflux ''beam-flux'' [photons/sec/mm{{{^}}}2]
  * beampol_factor ''polarisation-factor-of-beam'' [fraction]
  * beampol_direct ''direction-of-polarisation-with-respect-to-roation-axis'' [degrees]


wavelength 0.24344        # wavelength in Angstrom
beam_width 0.8            # Beam width (mm)
                          # If no beam width is specified it is assumed that the entire sample
                          # width is illuminated
beamflux  1e13            # Beam flux (Ph/s/mm2)
beampol_factor 1.0        # Beam polarisation factor: 1 = fully plane polarised, 0 = unpolarised
beampol_direct 0.0        # Direction of the normal to the plane of the primary beam polarisation
                          # with respect to the sample rotation axis (degrees) e.g. if the omega
                          # rotation axis is parallel to the laboratory z-axis the value is
                          # 0.0 degrees and if along y-axis it is 90.0 degrees
}}}

2.1.2 Detector specification 

Beam center on detector -
  * dety_center ''beam-center-in-y-direction'' [pixels]
  * detz_center ''beam-center-in-y-direction'' [pixels]

  E.g.
    dety_center 1023.5               # beamcenter, y in pixel coordinates
    detz_center 1023.5               # beamcenter, z in pixel coordinates


Detector pixel size -
  * y_size ''pixel-size-y-direction'' [mm]
  * z_size ''pixel-size-z-direction'' [mm]

  E.g.
    y_size         0.04677648        # Pixel size y (mm)
    z_size         0.04808150        # Pixel size z (mm)


Detector size -
  * dety_size ''detector-size-y-direction'' [pixels]
  * detz_size ''detector-size-z-direction'' [pixels]

  E.g.
    dety_size   2048.0               # detector y size (pixels)
    detz_size   2048.0               # detector z size (pixels)

Distance from sample to detector -
  * distance ''distance-sample-to-detector'' [mm]

  E.g.
    distance      55.0               # sample-detector distance (mm)

Detector tilts -
  * tilt_x ''tilt-of-detector-ccw-around-x-axis'' [radians]
  * tilt_y ''tilt-of-detector-ccw-around-y-axis'' [radians]
  * tilt_z ''tilt-of-detector-ccw-around-z-axis'' [radians]

  The order of the tilts is Rx*Ry*Rz, thus the tilt (rotation) 
  about z is performed first, then about y and finally about x.

  E.g.
    tilt_x 0.0
    tilt_y 0.01
    tilt_z 0.0 

OBS: If diffraction images are to be formed it is also possible to
     simulate detector point spread, background and noise. See below.

Detector orientation -
  Two things determine the detector orientation: 
    1. How it is mounted in the beam line setup  
    2. How the device server reads out the image.
  To get the image in the standard FABLE geometry a orientation matrix o 
  can be specifed:

    o = [o11 o21]
        [o21 o22]

  (see this document Geometry_version_1.0.7.pdf on the wiki)

  There are eight possible '''o''' matrices for the eight possible orientations.
  * o11 ''element-in-orientation-matrix'' [-1,0,1]
  * o12 ''element-in-orientation-matrix'' [-1,0,1]
  * o21 ''element-in-orientation-matrix'' [-1,0,1]
  * o22 ''element-in-orientation-matrix'' [-1,0,1]

  E.g.
    o11  1               # Frelon2k, Frelon4m detector orientation
    o12  0               #
    o21  0               # 
    o22 -1               #

  OBS: In principle the choice of matrix does not matter - 
       only if you want to get the images in the same orientation 
       as what you get on a specific setup.
       Remember the values given here has to be the same if 
       analysing the data in ImageD11.

Background -
  A constant background can be added to the diffraction images using:
  * bg ''number_of_counts''

  E.g.
    bg 100                # Add 100 counts to background

Noise -
  Random Poisson noise can be added to the diffraction images using:
  * noise ''flag'' [0= no noise, 1= add Poisson noise]

  E.g.
    noise 1               # Add Poisson noise


Detector point spread -
  A Gaussian detector point spread can be simulated in the diffraction
  images by specifying the FWHM of the point spread in pixels:
  * psf ''fwhm'' [in pixels]

  E.g.
    psf 2                 # Add Gaussian detector psf with a FWHM of 2 pixels

Spatial distortion -
  Spatial distortion of the detector can be taken into account 
  by specifying the corresponding spline. The spots will not be
  deformed due to the operation, but merely the center position
  of the peak is distorted.
  * spatial '''spline file'''

  E.g.
    spatial 'frelon4m.spline'          # Add spatial distortion 


2.1.3 Omega scan range, step size and speed 

Omega scan range -
  * omega_start ''start-omega-value-of-scan'' [degrees]
  * omega_end ''end-omega-value-of-scan'' [degrees]

  E.g.
    omega_start  -45.0             # Minimum Omega in range of interest (in deg)
    omega_end     45.0             # Maximum Omega in range of interest (in deg)

Omega step size -
  * omega_step ''omega_step_size_of_each_frame'' [degrees]

  E.g.
    omega_step     1.0             # Omega step size (in deg)


Omega rotation direction -
  * omega_sign ''omega_rotation_direction'' [+1/-1]

  E.g.
    omega_sign     1         # Sign of omega rotation (cw = +1, ccw = -1)

  OBS: The above omega specifications will create 
       (omega_end - omega_start)/omega_step images, and the first image 
       will be centered at omega_start + 0.5*omega_step. In the file header 
       Omega will be set to the midpoint of the rotation interval for the frame.

Wedge angle of omega axis - 
  The angle between the omega rotation axis and the z-laboratory axis in the
  plane of x and z, i.e. if 0 the rotation axis is perfectly aligned with the
  z-axis. Hence the rotation of the rotation axis about the y-axis 
  (left-handed).

  * wedge ''wedge_angle'' [degrees]

  E.g.
    wedge          0.023           # wedge of omega rotation axis (in deg)


2.2   CRYSTAL/GRAIN PARAMETERS

Grain/crystal number in sample to be simulated.

  * no_grains ''number-of-simulated-grains''
   The number is the total of all grains summed over all phases to be
   simulated. This number needs to match the number of e.g. 
   U_grains_X keywords

  E.g.
    no_grains 10   

Grain phase -
  The phase of the individual grains can be specified or appointed by
  PolyXSim. If you want to let the PolyXSim appoint which grain belongs
  to which phase the following keyword can be used.

  * gen_phase ''flag'' [0= do not, or 1= do] ''phase-id no-of-grains-of-this-id'' .... ''phase-id no-of-grains-of-this-id''
  
  So if choosing to randomly appoint the phase a list of phase id's each 
  succeded by the number of grains of the present phase. 
  E.g. (example with three phases with id'd 0,1 and 2)

    gen_phase 1 0 10 1 20 2 5 
    10 grains with phase 0, 
    20 grains having phase 1, and
    5 grains with phase 2

  Naturally the number of grains should match the number of grains specified
  with keyword 'no_grains'

  Alternatively the phase id of each phase can be specified
  * phase_grains_X ''phase-id-of-grain''
  where X is the grain id number.

Grain orientations - 
  Can be either randomly generated, or specific orientation matrices
  can be input by the user.

  * gen_U ''flag'' [0= do not, or 1= do]
  If  flag = 1 then 'no_grains' random orientations will be generated

  E.g.
    gen_U 1              # Generate orientations

  If  gen_U is set to 0 the orientations have to be provided by the user
  * U_grains_''X'' ''U11 U12 U13 U21 U22 U23 U31 U32 U33 
  X needs to be integer - its used to make certain that a grain orientation 
  is correctly matched with its position, size etc.

  E.g. 
    U_grains_0 -0.888246 0.411253 -0.204671 -0.201101 -0.748709 -0.631659 -0.413011 -0.519909 0.747741
    U_grains_1 -0.158282 -0.986955 0.029458 -0.929214 0.158978 0.333597 -0.333929 0.025430 -0.942255
    ..........


Grain positions -
  Can be either randomly generated or specific positions can be input by 
  the user.

  * gen_pos value1 [0= do not, or 1= do] value2 [0= all at (0,0,0), 1= generate randomly within box or cylinder]

  If ''flag1''=1 ''no_grains'' random positions will be generated

  E.g.
    gen_pos 1 1        # Generate random positions within box or cylinder


  OBS: The function of gen_pos is dependent on other keywords (or lack
       of keywords). For generation of a position different from (0,0,0)
        one of the keywords sample_cyl or sample_xyz should be given in 
        order define the borders of the sample area.

  * pos_grains_X ''x y z'' [mm]
  X needs to be integer - its used to make certain that a grain position
  is correctly matched with its orientation, size etc.''

  E.g.
    pos_grains_0 0 0 0
    pos_grains_1 0.01 -0.05 0.2
    .......

Sample shape and dimensions - 
  The sample can be specified to have either cylindrical or box shape:
  * sample_cyl ''diameter height'' (dimensions given in mm)
  or
  * sample_xyz x_dimension y_dimension z_dimension (all in mm)

  OBS: Only one of sample_cyl and sample_xyz can be given.
       If no beam_width is given it is assumed that the entire
       width of the sample is illuminated.

  E.g.
    sample_cyl 0.8 0.1             # Cylindrical sample shape
      or
    sample_xyz 0.5 0.5 0.5         # Box shaped sample


Grain strains -
  Can be either randomly generated, or specific strains can be input
  by the user. Note that the strain tensor is given in the Cartesian
  grain coordinate system, which for each grain is related to the
  overall sample system via the grain specific orienation matrix U.

  * gen_eps flag [0= do not, or 1= do] mean-value-for-diagonal-elemets-of-strain-tensor 
                                        spread-for-diagonal-elements-of-strain-tensor
                                        mean-value-for-offdiagonal-elemets-of-strain-tensor
                                        spread-for-offdiagonal-elements-of-strain-tensor
  If flag = 1 then no_grains strain tensors with with elements from a normal distribution with the specified mean and spread will be generated

  E.g.
    gen_eps 1 0 0.001 0 0         # Generate random diagonal strain tensors

    OBS: if a multiphase material is simulated using the above keyword the
         strain for all grain independt of phase will be generated with the
         same distribution. It is also possible to have different distributions
         for every phase.

  * gen_eps_phase_Y has the same entries as gen_eps given above
  Y being the phase number id

  E.g.
    gen_eps_phase_0 1 0 0.001 0 0    # Generate random diagonal strain tensors
    gen_eps_phase_1 1 0 0.02 0 0.01

  The strain tensors can also been specifically input for every grain

  * eps_grains_X eps11 eps12 eps13 eps22 eps23 eps33
  X needs to be integer - its used to make certain that the strain tensor 
  of the grain is correctly matched with its position, size etc.

  E.g.
    eps_grains_0 0.001 0.0015 -0.005 0 0 0
    eps_grains_1 0.001 -0.005 0.002 0.006 -0.005 -0.001
    .......

Grain sizes -
  Again these can either be user supplied or generated by PolyXSim.
  The grain sizes will be simulated having a log-normal distribution 
  with a specified median grain size and optionally the distribution
  tails can be cut off. If only one phase is to be simulated or one
  wishes to use the same grain distribution for all structural phase
  the following keyword can be used to specific the distribution

  * gen_size ''flag'' [0= do not, or 1= do] ''median-grain-size-of-distribution'' [mm] ''minimum-grain-size'' [mm] ''maximum-grain-size'' [mm]

  OBS: if value ''median-grain-size-of-distribution'' is negative
       the grain size of all grains will be the value of the median.

  E.g.
    gen_size 1 0.05 0.01 0.25 

  Different grain size distributions can be used for the different phase
  (if more than one are present). This is specified as follows

  * gen_size_phase_Y ''flag'' [0= do not, or 1= do] ''median-grain-size-of-distribution'' [mm] ''minimum-grain-size'' [mm] ''maximum-grain-size'' [mm]
  where Y again is the phase number id. [[br]]
  E.g.
    So the input can look like this for two or more phases,
  
    gen_size_phase_0 1 0.05 0.010 0.25   
    gen_size_phase_1 1 0.02 0.005 0.05   
    .....


  Or the grain size of each grain can be specified

  * size_grains_''X'' ''grain-diameter'' [mm]

  E.g.
    size_grains_0 0.04
    size_grains_1 0.06
    ......

Structural parameters -

  It is possible to simulate both mono- and multiphase polycrystalline 
  samples. If there is no interest in the actual peak intensities - only
  the unit cell and space group have to be specified.

  * unit_cell_phase_Y a [AA] b [AA] c [AA] alpha [deg] beta [deg] gamma [deg] 
    and
  * sgno_phase_Y ''number-of-space-group'', Y being the phase number id
    or
  * sgname_phase_Y ''name-of-space-group'', Y being the phase number id
  Y being the phase number id

  Presently only the standard space groups can be used, i.e. P 21/n is
  for example not a possibility.

  E.g.
    unit_cell_phase_0  8.531200 4.832100 10.125000 90.000000 92.031000 90.000000
    sgno_phase_0  4         # space group number
      or
    sgname_phase_0 'P21'    # remember to put quotation marks around the string
  OBS:  if more phases the next set will then have the keywords
        unit_cell_phase_1, sgno_phase_1 etc.

  OBS2: If monophase materails are simulated the old keyword (i.e. 
        without _phase_Y ) can still used.

  If 'real' intensities are to be calculated the structural parameters can
  either be supplied as a [http://www.iucr.org/iucr-top/cif/index.html cif]
  file or a [http://www.wwpdb.org/docs.html pdb] file.

  * structure_phase_Y '''structure-file-name'''
  The file can either be a .pdb or a .cif file
  Y being the phase number id

  E.g.
    structure_phase_0 'glycine.cif'

  Again if a monophase is simulated the old keyword - structure_file - can 
  be used instead.

  If a structure file has information about unit cell and/or space group, 
  these parameter will be chosen over parameters introduced with the keywords
  unit_cell and/or sgno.

File names and formats -

 Directory to save output from PolyXSim -
  * direc '''directory-name'''
  If the specified directory does not exist it will be created.

  E.g.
    direc 'simulation_glycine'

Name stem -
  The base of all out put files
  * stem '''name-stem'''
  i.e. image files will get the names name-stem_frame0001.edf etc.

  E.g.
    stem 'glycine'

File formats to output -
  There is a number of files which will be made by default. Whether one
  likes it or not. But the simulated reflections can be out put in 
  different file formats if requested. What can be chosen
  1. '.edf' or '.tif' - presently the two supported diffraction image formats
  2. '.flt' - a peak file definitions on wiki:"imaged11 - file formats" format
  3. '.gve' - defined on the wiki:"imaged11 - file formats" g-vector file.
  4. '.ubi' - grain orientations as inv(U*B)
  5. '.par' - the input parameters for PolyXSim written in the par format 
              of ImageD11

  The output files are specified as one after the other in the following manner
  * output '''format1' 'format2' ......''

  E.g. 
    output '.edf' '.flt' '.gve' '.ubi' '.par'  

Peak intensities -

The total peak intensity can be either be a
  1. constant value, or
  2. based on the '''structure factor''' of the reflections

  The following keyword is used to control this choice

  * intensity_const ''constant_intensity_value'' [counts]
  if the value is zero (0) the structure factor squared are used to
  calculate the intensities. Otherwise the value given will be the 
  total intensity.

  In the calculation of the intensity the effects of
    1. the Lorentz factor and/or
    2. the beam polarisation factor
  can also be taken into account using
  * lorentz_apply ''flag'' [0=do not,1=do]
  * beampol_apply ''flag'' [0=do not,1=do]


Peak shapes -

  Currently only three different peak shapes/profiles are available
  * peakshape ''type''

  Depending on the chosen type more parameters can be added after ''type''. 
  Type can be
    * 0 - spike peaks
    * 1 - Gaussian peak shape (isotropic), and Gaussian
          rocking curve (spread along omega)
    * 2 - peak shape calculated from a grain orientation distribution function
  Below these three types are documented in more detail.

  0. Spike peak - A square 2-by-2 spike peak
  * peakshape 0

  E.g.
    peakshape 0

  1. Gaussian peaks

  Spot with one FWHM in the y and z detector directions (in pixels)
  and another in omega (degress)
  * peakshape 1 ''spot-full-width-half-maximum''[pixels] ''spot-rocking-curve'' [degrees]

  E.g.
    peakshape 1  2  0.2      

  2. Orientation spread dependent peaks 
  Use the 'real' orientation spread of the crystal (mosaicity). 
  OBS: This will not smear intensity in the 2theta direction. 
  NB: Only one common ODF can presently be used for all grains
      irrespective of phase.

  E.g.
    peakshape 2     # Make peak spread from orientation distribution function

  There is two possibilities of providing/defining an orientation
  distribution function (ODF)
  * odf_type ''type-code'' [1,2, or 3]

  1. The simplest one is to give to ODF as a isotropic Gaussian mosaic spread
  * odf_type 1
  the mosaicity is then given by the keyword:
  * mosaicity ''mosaic-spread'' [degrees]

  Optionally a scale of the ODF grid can be given - by default it is given
  a value of half the anguler size of a pixel having the smallest angular
  dimension.
  * odf_scale ''grid_voxel_size'' [degrees]

  E.g.
    odf_type 1
    mosaicity 0.2        # The mosaic spread in degrees

  If no odf_type or mosaicity the values above will be used by default.

  2. The other is to give the ODF as voxelated grid defined in Rodrigues space. 
     The ODF has to read from a file [wiki:"PolyXSim - odf_file format" format]
  * odt_type 2
  * odf_file '''odf-data-file''' [file format]
 
  E.g.
    odf_type 2
    odf_file 'my_odf.odf'
    """
