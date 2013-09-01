from GtBurst import aplpy
import matplotlib.pyplot as plt
import pyfits, time
from GtBurst import dataHandling

class InteractiveFt1Display(object):
  def __init__(self,ft1file,skyimage,figure):
    self.skyimage             = skyimage
    ft1                       = pyfits.open(ft1file)
    events                    = ft1['EVENTS'].data
    self.empty                = False
    if(len(events)==0):
      print("No events in FT1 file %s" %(ft1file))
      self.empty              = True
      #raise RuntimeError("No events in FT1 file %s" %(ft1file))
    pass
    
    #Read in the different classes of events
    self.trigtime             = dataHandling.getTriggerTime(ft1file)
    time                      = events.field("TIME")
    energy                    = events.field("ENERGY")
    if(not self.empty):
      self.tmin                 = min(time)-self.trigtime
      self.tmax                 = max(time)-self.trigtime
      self.energyMin            = min(energy)
      self.energyMax            = max(energy)
    else:
      self.tmin               = float(ft1['EVENTS'].header['TSTART'])
      self.tmax               = float(ft1['EVENTS'].header['TSTOP'])
      self.energyMin          = 100
      self.energyMax          = 1e7
    pass
    ra                        = events.field("RA")
    dec                       = events.field("DEC")
    bitmask                   = map(lambda x:bin(x),events.field('EVENT_CLASS'))
    ids                       = events.field("EVENT_ID")
    theta                     = events.field("THETA")
    zenith                    = events.field("ZENITH_ANGLE")
    try:
      self.transient          = filter(lambda x:x[-1][-1]=='1' and 
                                                  x[-1][-3]=='0' and 
                                                  x[-1][-4]=='0' and 
                                                  x[-1][-5]=='0',
                                                  zip(time,energy,ra,dec,ids,theta,zenith,bitmask))
    except:
      #No events?
      self.transient          = []
    try:
      self.source               = filter(lambda x:x[-1][-3]=='1' and 
                                                  x[-1][-4]=='0' and 
                                                  x[-1][-5]=='0',
                                                  zip(time,energy,ra,dec,ids,theta,zenith,bitmask))
    except:
      self.source             = []
    try:
      self.clean                = filter(lambda x:x[-1][-4]=='1' and 
                                                  x[-1][-5]=='0',
                                                  zip(time,energy,ra,dec,ids,theta,zenith,bitmask))
    except:
      self.clean              = []
    try:
      self.ultraclean           = filter(lambda x:x[-1][-5]=='1',
                                                   zip(time,energy,ra,dec,ids,theta,zenith,bitmask))
    except:
      self.ultraclean         = []
    pass
    
    ntransient                = len(self.transient)
    nsource                   = len(self.source)
    nclean                    = len(self.clean)
    nultraclean               = len(self.ultraclean)
    print("\nTransient class events:                     %s" %(ntransient+nsource+nclean+nultraclean))
    print("Source class events:                        %s" %(nsource+nclean+nultraclean))
    print("Clean class events:                         %s" %(nclean+nultraclean))
    print("Ultra clean class events:                   %s\n" %(nultraclean))
    
    self.forcedRefresh        = False
    self.pickerID             = None
    self.oldxmin              = -1e9
    self.oldxmax              = 1e9
    self.oldymin              = -1e9
    self.oldymax              = 1e9
    self.user_ra              = None
    self.user_dec             = None
    
    self.figure               = figure
    self.figure.clear()
    self.displayImage()
    self.initEventDisplay()
    self.displayEvents()
    self.figure.canvas.draw()
    self.connectEvents()
    ft1.close()
  pass
  
  def unbind(self):
    #Clear all bindings in figures
    #print("UNBINDING")
    self.rangerTimer.stop()
    self.timer.stop()
    self.figure.canvas.mpl_disconnect(self.pickerID)
    self.figure.canvas.mpl_disconnect(self.clickerID)
  pass
  
  def displayImage(self):
    self.image               = aplpy.FITSFigure(self.skyimage,convention='calabretta',
                                                figure=self.figure,
                                                subplot=[0.1,0.15,0.40,0.7],
                                                label='sky image')
    
    imageFits                = pyfits.open(self.skyimage)
    img                      = imageFits[0].data
    
    # Apply grayscale mapping of image
    if(not self.empty):
      skm                      = self.image.show_colorscale(cmap='gist_heat',vmin=0.1,
                                                 vmax=max(img.flatten()),
                                                 stretch='log')
    else:
      skm                      = self.image.show_colorscale(cmap='gist_heat',vmin=0,
                                                 vmax=0.1)
    imageFits.close()
    
    # Modify the tick labels for precision and format
    self.image.tick_labels.set_xformat('ddd.dd')
    self.image.tick_labels.set_yformat('ddd.dd')
    
    # Display a grid and tweak the properties
    try:
      self.image.show_grid()
    except:
      #show_grid throw an expection if the grid was already there
      pass
    self.figure.canvas.draw()
  pass
  
  def initEventDisplay(self):
    #This must be called just once!
    self.eventDisplay         = self.figure.add_axes([0.60,0.15,0.35,0.7],label='event display')
  pass
  
  def inRegion(self,ra,dec,xmin,xmax,ymin,ymax):
    #Transform in pixel coordinates then check if ra,dec is contained
    #in the provided rectangular region
    tr                        = self.image._ax1._wcs.wcs_sky2pix
    x,y                       = tr(ra,dec,1)
    if((x>=xmin and x<=xmax) and
       (y>=ymin and y<=ymax)):
       #Contained
       return True
    else:
       return False
    pass
  pass
  
  def displayEvents(self,xmin=-1,xmax=1e9,ymin=-1,ymax=1e9):
    
    #Filter data
    
    transient                 = filter(lambda x:self.inRegion(x[2],x[3],xmin,xmax,ymin,ymax),
                                                self.transient)
    source                    = filter(lambda x:self.inRegion(x[2],x[3],xmin,xmax,ymin,ymax),
                                                self.source)
    clean                     = filter(lambda x:self.inRegion(x[2],x[3],xmin,xmax,ymin,ymax),
                                                self.clean)
    ultraclean                = filter(lambda x:self.inRegion(x[2],x[3],xmin,xmax,ymin,ymax),
                                                self.ultraclean)
    #Events display
    self.eventDisplay.cla()
    self.eventDisplay.scatter(map(lambda x:x[0]-self.trigtime,transient),
                 map(lambda x:x[1],transient),s=10,label='Transient',color='grey',
                 picker=0)
    self.eventDisplay.scatter(map(lambda x:x[0]-self.trigtime,source),
                 map(lambda x:x[1],source),s=30,label='Source',color='red',
                 picker=0)
    self.eventDisplay.scatter(map(lambda x:x[0]-self.trigtime,clean),
                 map(lambda x:x[1],clean),s=30,label='Clean',color='blue',
                 picker=0)
    self.eventDisplay.scatter(map(lambda x:x[0]-self.trigtime,ultraclean),
                 map(lambda x:x[1],ultraclean),s=30,label='Ultra Clean',color='cyan',
                 picker=0)
    try:
      self.eventDisplay.set_yscale('log')
      self.eventDisplay.set_ylim([self.energyMin*0.8,self.energyMax*1.2])
      self.eventDisplay.set_xlim([self.tmin-0.3*abs(self.tmin),self.tmax+0.3*abs(self.tmax)])  
    except:
      #no events to display, restore linear mode otherwise "figure.canvas.draw()"
      #will fail (!)
      self.eventDisplay.set_yscale('linear')
      pass
    self.eventDisplay.set_ylabel("Energy (MeV)",fontsize='small')
    self.eventDisplay.set_xlabel("Time since trigger (s)",fontsize='small')
    #Put the legend on top of the figure
    legend                    = self.eventDisplay.legend(scatterpoints=1,
                                            ncol=2,labelspacing=0.05,
                                            loc='upper center',
                                            bbox_to_anchor=(0.5,1.20),
                                            fancybox=True)
    legend.get_title().set_fontsize('6')
    self.forcedRefresh        = True
    self.figure.canvas.draw()
    
    #Destroy the callback with previous data, and create a new one with the new data
    if(self.pickerID!=None):
      self.figure.canvas.mpl_disconnect(self.pickerID)
    pass
    self.pickerID             = self.figure.canvas.mpl_connect('pick_event',
                                                                lambda x:self.on_events_plot_click(x,
                                                                                     transient,
                                                                                     source,
                                                                                     clean,
                                                                                     ultraclean))
  pass
  
  def connectEvents(self):
    self.clickerID            = self.figure.canvas.mpl_connect('button_press_event',self.on_click)
    self.timer                = self.figure.canvas.new_timer(interval=200)
    self.timer.add_callback(self.keep_synch, self.image._ax1)
    self.rangerTimer          = self.figure.canvas.new_timer(interval=100)
    self.rangerTimer.add_callback(lambda x:self.clearRanger(x,self.rangerTimer),self.figure)
    self.timer.start()
    self.rangerTimer.start()
  pass
  
  def waitClick(self):
    #Put the window in waiting mode, waiting for a click on the sky image
    self.locking              = True
    self.figure.canvas.start_event_loop(0)
  pass
  
  def on_click(self,event):
    if(event.inaxes==None):
      #Click outside any plot, do nothing
      return
    pass
    
    if(event.inaxes.get_label()!='event display'):
      #This is a click on the sky image, get the corresponding ra,dec
      ax                      = self.image._ax1
      if(ax.get_navigate_mode()==None):
        ra,dec                = ax._wcs.wcs_pix2sky(event.xdata,event.ydata,1)
        self.user_ra          = ra[0]
        self.user_dec         = dec[0]
        self.figure.canvas.stop_event_loop()
      else:
        #We are in zoom mode, do nothing
        return
    else:
      #the click was not on the sky image. Do nothing
      return
    pass
  pass
  
  def clearRanger(self,event,rangerTimer):
    #Verify if the figure has been cleared. If so, remove all bindings
    if(len(self.figure.get_axes())==0):
      rangerTimer.stop()
      self.unbind()
    pass
  pass
  
  def keep_synch(self,event=None):
    ax                     = self.image._ax1
    
    nx                     = ax.get_xlim()
    ny                     = ax.get_ylim()
    ras, decs              = ax._wcs.wcs_pix2sky(nx,ny,1)
    xmin                   = nx[0]
    xmax                   = nx[1]
    ymin                   = ny[0]
    ymax                   = ny[1]
    if(xmin!=self.oldxmin or ymin!=self.oldymin or
       xmax!=self.oldxmax or ymax!=self.oldymax):
      self.displayEvents(xmin,xmax,ymin,ymax)
      self.oldxmin           = xmin
      self.oldxmax           = xmax
      self.oldymin           = ymin
      self.oldymax           = ymax
      #refresh, avoiding an infinite loop
      self.figure.canvas.draw()
    else:
      #Not changed
      return
  pass
  
  def on_events_plot_click(self,event,transient,source,clean,ultraclean):
    ind                 = event.ind[0]
    data_class          = event.artist.get_label().lower()
    
    if(data_class=='transient'):
      v                 = transient
    elif(data_class=='source'):
      v                 = source
    elif(data_class=='clean'):
      v                 = clean
    elif(data_class=='ultra clean'):
      v                 = ultraclean
    pass
    
    try:
      ras             = v[ind][2]
      decs            = v[ind][3]
      event_id        = v[ind][4]
      theta           = v[ind][5]
      zenith          = v[ind][6]
    except:
      print("Could not get Ra,Dec of your event. Please retry...")
      pass
    
    #Get the width and height (in deg) of the image display
    ax                     = self.image._ax1
    nx                     = ax.get_xlim()
    ny                     = ax.get_ylim()
    rass, decss            = ax._wcs.wcs_pix2sky(nx,ny,1)
    img_width              = min(abs(rass[0]-rass[1]),abs(decss[0]-decss[1]))
    radius_length          = 0.3
    
    self.forcedRefresh     = True
    self.image.show_circles([ras],[decs],[radius_length],
                           facecolor='white',layer="circle",alpha=0.8)
    self.forcedRefresh     = True
    self.image.add_label(ras-radius_length,decs-radius_length,
                         'id = %s\ntheta = %3.1f\nzenith = %3.1f' % (event_id,theta,zenith),
                         color='white',layer='event_id',
                         verticalalignment='top',size='small')

    self.forcedRefresh = True
    self.figure.canvas.draw()
  pass
pass
