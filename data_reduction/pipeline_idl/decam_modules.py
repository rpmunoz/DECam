class Worker_sex(multiprocessing.Process):

  def __init__(self, work_queue, result_queue):

    # base class initialization
    multiprocessing.Process.__init__(self)

    # job management stuff
    self.work_queue = work_queue
    self.result_queue = result_queue
    self.kill_received = False

  def run(self):
    while not self.kill_received:

      # get a task
      try:
        im_file, cat_file, seg_im_file = self.work_queue.get_nowait()
      except Queue.Empty:
        break

      fwhm=1.
      pixel_scale=0.263
      weight_type='MAP_WEIGHT'
      checkimage_type='NONE'
      checkimage_file='NONE'
      satur_level=4.3e5
      analysis_thresh=2.0
      detect_minarea=3
      detect_thresh=1.4
      phot_apertures=",".join(["%.2f" % x for x in 2*np.array((0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.))*fwhm/pixel_scale])
      filter_name='sex_config/gauss_3.0_7x7.conv'
      xml_name='ngvs_sex.xml'

      # the actual processing
      log_file="ngvsir_completeness_sex_thread%d.log" % i_thread
      for i in range(i_range[0],i_range[1]):
        command = "sex %s -c sex_config/cfht_wircam.sex -PARAMETERS_NAME sex_config/cfht_wircam_psfmodel.param -CATALOG_TYPE FITS_LDAC -CATALOG_NAME %s -SEEING_FWHM %.2f -WEIGHT_TYPE %s -WEIGHT_THRESH 0. -WEIGHT_IMAGE %s -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s -SATUR_LEVEL %d -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 30 -BACK_SIZE 250 -BACK_FILTERSIZE 3 -MASK_TYPE CORRECT -ANALYSIS_THRESH %.2f -DETECT_MINAREA %d -DETECT_THRESH %.2f -DEBLEND_MINCONT 0.0000001 -INTERP_TYPE ALL -INTERP_MAXXLAG 1 -INTERP_MAXYLAG 1 -FLAG_TYPE OR -FLAG_IMAGE %s -PHOT_AUTOPARAMS 2.3,4.0 -PHOT_FLUXFRAC 0.5 -PHOT_APERTURES %s -PIXEL_SCALE %.4f -FILTER Y -FILTER_NAME %s -WRITE_XML Y -XML_NAME %s -PSF_NAME %s -PSF_NMAX 1" % (mock_im_file[i], sex_cat_file[i], fwhm, weight_type, weight_file[i], checkimage_type, checkimage_file, satur_level, analysis_thresh, detect_minarea, detect_thresh, flag_file[i], phot_apertures, pixel_scale, filter_name, xml_name, psf_file[i] )
        print command
        with open(log_file, "a") as log:
          result=subprocess.call(command, stderr=log, stdout=log, shell=True)
#       log=subprocess.Popen(command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

        log.close()
        print "SExtractor thread: %d - iteration: %d is done!" % (i_thread, i)

      self.result_queue.put(id)



