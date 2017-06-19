#!/usr/bin/env python
# by Ran Novitsky Nof (ran.nof@gmail.com), 2015 @ BSL

# ***********************************************************************************
# *    Copyright (C) by Ran Novitsky Nof                                            *
# *                                                                                 *
# *    E2ReviewTool.py is free software: you can redistribute it and/or modify      *
# *    it under the terms of the GNU Lesser General Public License as published by  *
# *    the Free Software Foundation, either version 3 of the License, or            *
# *    (at your option) any later version.                                          *
# *                                                                                 *
# *    This program is distributed in the hope that it will be useful,              *
# *    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
# *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
# *    GNU Lesser General Public License for more details.                          *
# *                                                                                 *
# *    You should have received a copy of the GNU Lesser General Public License     *
# *    along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
# ***********************************************************************************

import argparse,sys,os,re,math,time
import matplotlib as mpl
mpl.use('QT4Agg')
from obspy import UTCDateTime
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import uic
from matplotlib.backend_bases import NavigationToolbar2, Event
from matplotlib.backends.backend_qt4agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from matplotlib import cm
from numpy import array,dtype,int,str,float,loadtxt,append,log10,where
import pyproj
geo = pyproj.Geod(ellps='WGS84')
# util class for Open Street Map
from osm import OSM as osm

parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description='''E2ReviewTool - Review ElarmS E2 log file''',
         epilog='''Created by Ran Novitsky Nof (ran.nof@gmail.com), 2015 @ BSL''')
parser.add_argument('-i',metavar='file',nargs='+',help='input E2 log file(s)',type=str)
parser.add_argument('-r',metavar='file',nargs='+',help='''input reference csv file
csv - comma separated value in GII DB format [epiid,ml,typ,lat,lon,abs_ot_d,abs_ot_t,refDI].
epiid: Event Catalog ID
ml: Local Magnitude
typ: Event type code
lat: Latitude
lon: Longitude
abs_ot_d: Origin time Date
abs_ot_t: Origin time time
refDI: E2 reference ID (or None)
 ''',type=str)
parser.add_argument('-b',metavar='bounds',nargs=4,help='Region bounding box (west east south north)',type=float,default=(-180.0,180.0,-90,90))

FONTSIZE=8
VERBOSE=False # printout message
GRIDON=False # grid on or off [True | False]
OSMTILEURL="http://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/" # where to read map tiles from
OSMTILEARCHIVE='tiles1' # where to save downloaded tiles
OSMTILEPAT = "{Z}/{Y}/{X}.png"
TILEARCHIVE='../ElViS/tiles' # local map tiles

##################### Some matplotlib Black Magic ##########################################################
## This part will redirect some matplotlib toolbar and canvas callbacks
## adjusting navigation bar for capturing after zoom/pan events

# reference to original matplotlib toolbar functions
canvasHome = NavigationToolbar2.home
canvasBack = NavigationToolbar2.back
canvasForward = NavigationToolbar2.forward
canvasRelease_pan = NavigationToolbar2.release_pan
canvasRelease_zoom = NavigationToolbar2.release_zoom

# new toolbar functions with a callback signal (defined by a)
# each function will call the original toolbar function and emite a callback
# signals can be redirected to a new function at a later stage using mpl_connect.
def new_home(self, *args, **kwargs):
  a = 'after_home_event' # this is the signal called after processing the event
  event = Event(a, self)
  canvasHome(self, *args, **kwargs) # call original matplotlib toolbar function
  self.canvas.callbacks.process(a, event) # process the signal
def new_back(self, *args, **kwargs):
  a = 'after_back_event'
  event = Event(a, self)
  canvasBack(self, *args, **kwargs)
  self.canvas.callbacks.process(a, event)
def new_forward(self, *args, **kwargs):
  a = 'after_forward_event'
  event = Event(a, self)
  canvasForward(self, *args, **kwargs)
  self.canvas.callbacks.process(a, event)
def new_release_pan(self, evt):
  a = 'after_release_pan_event'
  event = Event(a, self)
  canvasRelease_pan(self,evt)
  self.canvas.callbacks.process(a, event)
def new_release_zoom(self, evt):
  a = 'after_release_zoom_event'
  event = Event(a, self)
  canvasRelease_zoom(self,evt)
  self.canvas.callbacks.process(a, event)
# change toolbar functions to the new functions
NavigationToolbar2.home = new_home
NavigationToolbar2.back = new_back
NavigationToolbar2.forward = new_forward
NavigationToolbar2.release_pan = new_release_pan
NavigationToolbar2.release_zoom = new_release_zoom
######################### End of matplotlib black magic ########################################

def E2pdmag(logpd,distkm):
  'Calculate magnitude as in ElarmS (Kuyuk and Allen, 2013)'
  return 5.39+1.23*logpd+1.38*log10(distkm)

def Sadehpdmag(logpd,distkm):
  'Calculate magnitude as in Sadeh, Ziv, Wust-Bloch, 2014 GJI'
  return 5.7935+(logpd+log10(distkm))/1.041

def concatenateCSV(CSVs):
  giidbtype = dtype([('epiid','|S50'),('ml',float),('typ',int),('lat',float),('lon',float),('abs_ot_d','|S50'),('abs_ot_t','|S50')])
  rettype = dtype([('EID','|S50'),('ot','|S50'),('lat',float),('lon',float),('depth',float),('mag',float),('refID','|S50'),('type',int)])
  ret = array([],dtype=rettype)
  for f in CSVs:
    try:
      fdata = loadtxt(f,delimiter=',',dtype=giidbtype,skiprows=1)
    except ValueError:
      raise ValueError('%s is not a valid scv file'%(f))
    ret = append(ret,array([(v[0],'T'.join([v[5],v[6]])+'Z',v[3],v[4],10.0,v[1],None,v[2]) for v in fdata],dtype=ret.dtype))
  return dict(zip(ret['EID'],ret))

def concatenateNewLOG(LOGs):
  originDtype = dtype([('EID','|S50'),('ver',int),('ot','|S50'),('lat',float),('lon',float),('depth',float),('mag',float),('nT',int),('nS',int),('Spercent',float),('alert',int),('alertime','|S50'),('first',int)])
  triggerDtype = dtype([('EID','|S50'),('ver',int),('ID','|S50'),('net','|S50'),('sta','|S50'),('loc','|S50'),('chn','|S50'),('trigT','|S50'),('lat',float),('lon',float),('pdmag',float),('logpd',float),('pdSNR',float),('updm',int),('distkm',float),('azimuth',float),('TTErr',float)])
  origins = array([],originDtype)
  triggers = array([],triggerDtype)
  originsdict = {}
  triggerdict = {}
  for f in LOGs:
    try:
      with open(f,'r') as log:
        for line in log:
          if re.match(".+E:I:[ F].+",line):
            try: # for older logs than July 22, 2015
              timeStamp,Eid,ver,lat,lon,depth,mag,otTxt,latu,lonu,depu,magu,timeu,lk,nTb,nSb,nT,nS,ave,rms,fitok,splitok,near,statrig,active,inact,nsta,percnt,prcntok,mindist,maxdist,distok,azspan,Mok,nSok,Lok,Tdif,Tok,Aok,Ast,atimeTxt = line.strip().split()
            except ValueError: # for logs since July 22, 2015
              timeStamp,Eid,ver,lat,lon,depth,mag,otTxt,latu,lonu,depu,magu,timeu,lk,nTb,nSb,nT,nS,ave,rms,fitok,splitok,near,statrig,active,inact,nsta,percnt,prcntok,mindist,maxdist,distok,azspan,Mok,nSok,Lok,Tdif,tpave,pdave,Tok,Azok,Aok,Ast,atimeTxt = line.strip().split()
            if not int(Eid)>0: continue
            origin = array([(Eid,ver,otTxt,lat,lon,depth,mag,nT,nS,percnt,Ast,atimeTxt,0)],dtype=originDtype)
            if re.match(".+E:I:F.+",line): origin['first'] = 1
            origins = append(origins,origin)
          if re.match(".+E:I:T: .+",line):
            timeStamp,Eid,ver,update,order,sta,chn,net,loc,lat,lon,trigger_time,log_taup,taup_snr,log_pd,pd_snr,log_pv,pv_snr,pa,pa_snr,assoc,tpmag,utpm,pdmag,updm,uch,ukm,upd,ups,utp,uts,distkm,azimuth,tterr = line.strip().split()[:34]# parse line
            if not int(Eid)>0: continue
            loc = loc.replace('--','')
            triggers = append(triggers,array([(Eid,ver,'.'.join([net,sta,loc,chn]),net,sta,loc,chn,trigger_time,lat,lon,pdmag,log_pd,pd_snr,updm,distkm,azimuth,tterr)],dtype=triggerDtype))
    except Exception as E:
      sys.exit('Error in %s:\n%s\n in line:\n%s'%(f,str(E),line))
  origins.sort(order=['EID','ver'])
  for origin in origins:
    if not origin['EID'] in originsdict: originsdict[origin['EID']] = []
    originsdict[origin['EID']].append(origin)
  for trigger in triggers:
    ID = '-'.join([trigger['EID'],str(trigger['ver'])])
    if not ID in triggerdict: triggerdict[ID]=[]
    triggerdict[ID].append(trigger)
  return originsdict,triggerdict

def concatenateLOG(LOGs):
  originDtype = dtype([('EID','|S50'),('ver',int),('ot','|S50'),('lat',float),('lon',float),('depth',float),('mag',float),('nT',int),('nS',int),('Spercent',float),('alert',int),('alertime','|S50'),('first',int),('OK',int)])
  triggerDtype = dtype([('EID','|S50'),('ver',int),('ID','|S50'),('net','|S50'),('sta','|S50'),('loc','|S50'),('chn','|S50'),('trigT','|S50'),('lat',float),('lon',float),('pdmag',float),('logpd',float),('pdSNR',float),('updm',int),('tpmag',float),('logtp',float),('tpSNR',float),('utpm',int),('distkm',float),('azimuth',float),('TTErr',float)])
  origins = array([],originDtype)
  triggers = array([],triggerDtype)
  originsdict = {}
  triggerdict = {}
  line = '--- None ---'
  for f in LOGs:
    try:
      with open(f,'r') as log:
        for line in log:
          line = re.sub('(....)\/(..)\/(..) (..:)','\\1-\\2-\\3T\\4',line)# adjust time stamps
          if re.match(".+E:I:[ F].+",line):
            try:# for older logs than July 22, 2015
              timeStamp,Eid,ver,lat,lon,depth,mag,otTxt,latu,lonu,depu,magu,timeu,lk,nTb,nSb,nT,nS,ave,rms,fitok,splitok,near,statrig,active,inact,nsta,percnt,prcntok,mindist,maxdist,distok,azspan,Mok,nSok,Lok,Tdif,Tok,Aok,Ast,atimeTxt = line.strip().split()
            except ValueError: # for logs since July 22, 2015
              timeStamp,Eid,ver,lat,lon,depth,mag,otTxt,latu,lonu,depu,magu,timeu,lk,nTb,nSb,nT,nS,ave,rms,fitok,splitok,near,statrig,active,inact,nsta,percnt,prcntok,mindist,maxdist,distok,azspan,Mok,nSok,Lok,Tdif,tpave,pdave,Tok,Azok,Aok,Ast,atimeTxt = line.strip().split()
            if not int(Eid)>0: continue
            OK = int(Mok)<<3 | int(nSok)<<2 | int(Lok)<<1 | int(Tok)
            origin = array([(Eid,ver,otTxt,lat,lon,depth,mag,nT,nS,percnt,Ast,atimeTxt,0,OK)],dtype=originDtype)
            if re.match(".+E:I:F.+",line): origin['first'] = 1
            origins = append(origins,origin)
          if re.match(".+E:I:T: .+",line):
            try: # for older logs than July 22, 2015
              timeStamp,Eid,ver,update,order,sta,chn,net,loc,lat,lon,trigger_time,log_taup,taup_snr,log_pd,pd_snr,log_pv,pv_snr,pa,pa_snr,assoc,tpmag,utpm,pdmag,updm,uch,ukm,upd,ups,utp,uts,distkm,azimuth,tterr = line.strip().split()[:34]# parse line
              if not '.' in log_taup: raise ValueError
            except ValueError: # for logs since July 22, 2015
              timeStamp,Eid,ver,update,order,sta,chn,net,loc,lat,lon,trigger_time,rsmp,tsmp,log_taup,taup_snr,dsmp,log_pd,pd_snr,assoc,tpmag,utpm,pdmag,updm,uch,ukm,upd,ups,utp,uts,tel,tsec,distkm,azimuth,tterr,plen,sps,toffset,arrtime,protime,fndtime,quetime,sndtime,e2time,buftime,alert = line.strip().split() # parse line
            if not int(Eid)>0: continue
            loc = loc.replace('--','')
            if log_taup=='NA': log_taup=-9999
            if taup_snr=='NA': taup_snr=-9999
            trigger_time= trigger_time.replace('/','-')
            triggers = append(triggers,array([(Eid,ver,'.'.join([net,sta,loc,chn]),net,sta,loc,chn,trigger_time,lat,lon,pdmag,log_pd,pd_snr,updm,tpmag,log_taup,taup_snr,utpm,distkm,azimuth,tterr)],dtype=triggerDtype))
          if re.match(".+\|H: .+",line):
            timeStamp,Eid,ver,atimeTxt,otTxt,mag,foo,foo,foo,foo,lat,lon,foo,nT,nS,percnt,foo,foo,foo,oterr,foo,Ast = line.strip().split()
            origin = array([(Eid,ver,otTxt,lat,lon,8.0,mag,nT,nS,percnt,(Ast=='yes'),atimeTxt,0,15)],dtype=originDtype)
            if not len(origins['first'][where(origins['EID']==Eid)])>0 and origin['alert']: origin['first'] = 1
            origins = append(origins,origin)
            if int(ver)>0:
              trigUpdate = triggers[where(triggers['EID']==Eid)]
              trigUpdate = trigUpdate[where(trigUpdate['ver']==int(ver)-1)]
              for trigup in trigUpdate:
                trig = trigup.copy()
                trig['ver']=int(ver)
                triggers = append(triggers,trig)
            # update distance,azimuth
            indx  =  where(triggers['EID']==Eid)[0]
            indx = indx[where(triggers[indx]['ver']==int(ver))]
            for i in indx:
              azimuth,dist = geo.inv(float(lon),float(lat),triggers[i]['lon'],triggers[i]['lat'])[-2:]
              triggers[i]['distkm'] = dist/1000.0
              triggers[i]['azimuth'] = azimuth
          if re.match(".+\|L[0-9]+: .+",line):
            timeStamp,Eid,sta,net,chn,loc,lat,lon,trigger_time,log_taup,taup_snr,log_pd,pd_snr,log_pv,pv_snr,pa,pa_snr,assoc,tpmag,utpm,pdmag,updm,distkm = line.strip().split()
            ver = int(timeStamp.split('L')[1][:-1])
            loc = loc.replace('--','')
            distkm,azimuth,tterr = -9999,-9999,-9999
            if taup_snr=='NA': taup_snr=-9999
            if log_taup=='NA': log_taup=-9999
            triggers = append(triggers,array([(Eid,ver,'.'.join([net,sta,loc,chn]),net,sta,loc,chn,trigger_time,lat,lon,pdmag,log_pd,pd_snr,updm,tpmag,log_taup,taup_snr,utpm,distkm,azimuth,tterr)],dtype=triggerDtype))
    except Exception as E:
      sys.exit('Error in %s:\n%s\n in line:\n%s'%(f,str(E),line))
  origins.sort(order=['EID','ver'])
  for origin in origins:
    if not origin['EID'] in originsdict: originsdict[origin['EID']] = []
    originsdict[origin['EID']].append(origin)
  for trigger in triggers:
    ID = '-'.join([trigger['EID'],str(trigger['ver'])])
    if not ID in triggerdict: triggerdict[ID]=[]
    triggerdict[ID].append(trigger)
  return originsdict,triggerdict


class CatalogEvent(QTreeWidgetItem):
  def __init__(self,EID='',ot='',lat=0,lon=0,depth=0,mag=0,refID='',typ=-1):
    QTreeWidgetItem.__init__(self,[EID,str(ot),'%.1f'%mag,'%.3f'%lat,'%.3f'%lon]+['']*6+[refID,'%d'%typ])
    self.EID = EID
    if not ot: ot=0
    if not lat:
      lat=0
      self.setText(3,'')
    if not lon:
      lon=0
      self.setText(4,'')
    if not depth: depth=0
    if not mag:
      mag=0
      self.setText(2,'')
    self.ot = UTCDateTime(ot)
    self.lat = float(lat)
    self.lon = float(lon)
    self.depth = float(depth)
    self.mag = float(mag)
    self.typ =typ
    self.refID = refID
    self.origins = []
    self.preffered = None
    if lat and lon:
      self.l = mpl.lines.Line2D([self.lon],[self.lat],color='k',marker='o',label=' '.join([EID,str(ot)]))
    else:
      self.l = mpl.lines.Line2D([],[],color='k',marker='o',label=' '.join([EID,str(ot)]))
    self.l.item = self
  def addOrigin(self,origin):
    def verSort(x,y):
      if x.ver==y.ver: return 0
      if x.ver>y.ver:
        return 1
      else:
        return -1
    origin.setEvent(self)
    self.origins.append(origin)
    km = geo.inv(origin.lon,origin.lat,self.lon,self.lat)[-1]/1000.
    Secs = origin.ot-self.ot
    if origin.first:
      self.setText(5,str(origin.ot))
      self.setText(6,'%.1f'%(origin.mag))
      self.setText(7,'%.2f'%(origin.lat))
      self.setText(8,'%.2f'%(origin.lon))
      if self.lat and self.lon:
        self.setText(9,'%.2f'%(km))
        self.setText(10,'%.2f'%(Secs))
      self.setText(11,str(origin.EID))
      self.preffered = origin
    self.origins.sort(verSort)


class Origin(QTreeWidgetItem):
  def __init__(self,EID,ver,ot,lat,lon,depth=10,mag=0,nT=0,nS=0,percnt=0,Aok=0,AT=0,first=0,OK=15):
    QTreeWidgetItem.__init__(self,['-'.join([EID,'%d'%ver]),str(AT),str(ot),'%.1f'%mag,'%.3f'%lat,'%.3f'%lon,'','','%d'%nT,'%d'%nS,'%.2f'%percnt,['no','yes'][Aok],'%d'%OK])
    if int(first):
      i = QIcon.fromTheme('emblem-default',QIcon(":/emblem-default.png"))
      self.setIcon(0,i)
    self.EID = EID
    self.ver = ver
    self.ot = UTCDateTime(ot)
    self.lat = float(lat)
    self.lon = float(lon)
    self.depth = float(depth)
    self.mag = float(mag)
    self.km = 0
    self.Secs = 0
    self.nT = int(nT)
    self.nS = int(nS)
    self.percnt = float(percnt)
    self.Aok = Aok
    self.OK = OK
    self.AT = UTCDateTime(AT)
    self.first=int(first)
    self.l = mpl.lines.Line2D([self.lon],[self.lat],color='k',marker='o',label=' '.join([EID,str(ver),str(ot)]))
    self.l.item = self
    self.triggers = []
    self.event = None
  def addTrigger(self,trig):
    self.triggers.append(trig)
    trig.origin = self
  def addTriggers(self,trigs):
    [self.triggers.append(trig) for trig in trigs]
  def setEvent(self,event):
    self.event = event
    if event.lat and event.lon:
      self.km = geo.inv(event.lon,event.lat,self.lon,self.lat)[-1]/1000.
      self.setText(6,"%.2f"%self.km)
      self.Secs = self.ot-event.ot
      self.setText(7,"%.2f"%self.Secs)

class StationMagnitude(QTreeWidgetItem):
  def __init__(self,EID,ver,ID,net,sta,loc,chn,tt,lat,lon,pdmag,logpd=0,pdsnr=0,updm=0,tpmag=0,logtp=0,tpsnr=0,utpm=0,dist=0,azim=0,tterr=0):
    QTreeWidgetItem.__init__(self,['-'.join([EID,str(ver)]),ID,str(tt),'%.3f'%lat,'%.3f'%lon,'%.1f'%dist,'%.2f'%tterr,'%.1f'%pdmag,['n','y'][updm],'%.2f'%logpd,'%.1f'%pdsnr,'%.1f'%tpmag,['n','y'][utpm],'%.2f'%logtp,'%.1f'%tpsnr,'%d'%azim])
    self.EID = EID
    self.ver = ver
    try:
      self.tt = UTCDateTime(tt)
    except:
      print >> sys.stderr,'Error with trigger timetamp:',EID,ver,tt
      sys.exit()
    self.lat = float(lat)
    self.lon = float(lon)
    self.logpd = float(logpd)
    self.pdsnr = float(pdsnr)
    self.pdmag = float(pdmag)
    self.updm = float(updm)
    self.logtp = float(logtp)
    self.tpsnr = float(tpsnr)
    self.tpmag = float(tpmag)
    self.utpm = float(utpm)
    self.dist = float(dist)
    self.azim = float(azim)
    self.tterr = float(tterr)
    self.ID = ID
    self.l = mpl.lines.Line2D([self.lon],[self.lat],color='k',marker='^',label=self.ID)
    self.l.item = self
    self.origin = None

class ConnectedLineEdit(QLineEdit):
  'QLineEdit class for connected limits (like max/min magnitude)'
  def __init__(self,val,minval,maxval,deci=2):
    QLineEdit.__init__(self)
    self.val = float(val)
    self.setText(val)
    self.minval=minval
    self.maxval=maxval
    self.validator1 = QDoubleValidator()
    self.validator1.setRange(minval,maxval,deci)
    #self.setValidator(validator)
    self.editingFinished.connect(self.validate)
    self.setToolTip('%0.2lf <= Value >= %0.2lf'%(minval,maxval))
    self.connectedLineEdit = None
    self.connectedLineEditVal = None
    self.connectedVal = None
    self.limDirection = None
    self.connected = False
  def setConnected(self,widget,param):
    assert param in ['top','bottom']
    self.connectedLineEdit = widget
    self.connectedLineEditParam = param # should be top or bottom
    self.connected=True
  def validate(self):
    if not self.validator1.validate(self.text(),2)==(2,2): # reset value if not valid
      self.setText(str(self.val))
      return 0
    if self.connected: # update connected widget limit or reset
      param = self.connectedLineEditParam
      widget = self.connectedLineEdit
      if param=='top':
        o = '<='
      else:
        o = '>='
      if eval(str(widget.text())+o+str(self.text())):
        self.val = float(self.text())
        eval('widget.validator1.set'+param.capitalize()+'('+str(self.val)+')')
        minval,maxval = (widget.validator1.bottom(),widget.validator1.top())
        widget.setToolTip("%0.2lf <= Value >= %0.2lf"%(minval,maxval))
        self.emit(SIGNAL('ok'))
        return 1
      else:
        self.setText(str(self.val))
        return 0
    self.val = float(self.text())
    self.emit(SIGNAL('ok'))
    return 1

# map area dialog
class bboxForm(QDialog):
  def __init__(self,parent=None):
    QDialog.__init__(self,parent=parent)
    self.setWindowTitle('SCxmlDiff - Set area rectangle')
    vbox = QVBoxLayout(self)
    w = QWidget()
    grid = QGridLayout(w)
    vbox.addWidget(w)
    westlabel = QLabel('West')
    W = ConnectedLineEdit('-180',-180,180)
    grid.addWidget(westlabel,2,1)
    grid.addWidget(W,2,2)
    eastlabel = QLabel('East')
    E = ConnectedLineEdit('180',-180,180)
    grid.addWidget(eastlabel,2,5)
    grid.addWidget(E,2,6)
    northlabel = QLabel('North')
    N = ConnectedLineEdit('90',-90,90)
    grid.addWidget(northlabel,1,3)
    grid.addWidget(N,1,4)
    southlabel = QLabel('South')
    S = ConnectedLineEdit('-90',-90,90)
    grid.addWidget(southlabel,3,3)
    grid.addWidget(S,3,4)
    self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel,Qt.Horizontal, self)
    vbox.addWidget(self.buttons)
    self.buttons.accepted.connect(self.accept)
    self.buttons.rejected.connect(self.reject)
    self.W = W
    self.E = E
    self.S = S
    self.N = N
  def setLims(self,w,e,s,n):
    self.W.setText(str(w))
    self.E.setText(str(e))
    self.S.setText(str(s))
    self.N.setText(str(n))
  def getLims(self):
    w = self.W.text().toDouble()[0]
    e = self.E.text().toDouble()[0]
    s = self.S.text().toDouble()[0]
    n = self.N.text().toDouble()[0]
    return w,e,s,n
  def validate(self):
    w,e,s,n = self.getLims()
    if w>=e and s>=n:
      QMessageBox.warning(self,'E2ReviewTool - Error','West value should be lower than East value.\nSouth value should be lower than North value.')
      return 0
    elif w>=e:
      QMessageBox.warning(self,'E2ReviewTool - Error','West value should be lower than East value.')
      return 0
    elif s>=n:
      QMessageBox.warning(self,'E2ReviewTool - Error','South value should be lower than North value.')
      return 0
    else:
      return 1

class Overlay(QWidget):
  'An overlay layer when updating'
  def __init__(self, parent = None):
    QWidget.__init__(self, parent)
    # Make transparent
    palette = QPalette(self.palette())
    palette.setColor(palette.Background, Qt.transparent)
    self.setPalette(palette)
    self._active = False
    self.timer=None
  def paintEvent(self, event):
    painter = QPainter()
    painter.begin(self)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.fillRect(event.rect(), QBrush(QColor(255, 255, 255, 127)))
    painter.setPen(QPen(Qt.NoPen))
    for i in range(6):
      if (self.counter / 5) % 6 == i:
        painter.setBrush(QBrush(QColor(127 + (self.counter % 5)*32, 127, 127)))
      else:
        painter.setBrush(QBrush(QColor(127, 127, 127)))
      painter.drawEllipse(
      self.width()/2 + 30 * math.cos(2 * math.pi * i / 6.0) - 10,
      self.height()/2 + 30 * math.sin(2 * math.pi * i / 6.0) - 10,20, 20)
    painter.end()
  def showEvent(self, event):
    if self._active: return
    self._active=True
    self.timer = self.startTimer(50)
    self.counter = 0
  def timerEvent(self, event):
    self.counter += 1
    #if self.counter == 60: self.counter=0
    self.update()
  def hideEvent(self,event):
    if not self._active: return
    self.killTimer(self.timer)
    self._active=False

class AppForm(QMainWindow):
  'QT form for application'
  def __init__(self, splash,args,parent=None):
    self.starttime=UTCDateTime(0) # start time of data
    self.endtime=UTCDateTime() # end time of data
    self.bbox=args.b # geographical bounding box [w,e,s,n]
    self.preferredOrigin='p' # prefferred should be 'p': preferred origin ; 'f': first origin ; 'l': last origin
    self.reportedOnly=False # use only reportded events
    self.minmag=0 # minimum magnitude
    self.maxmag=10.0 # maximum magnitude
    self.deltaT=60.0 # absolute time difference in seconds to associate events
    self.deltaR=300.0 # absolute distance difference in km to associate events
    self.deltaM=None # absolute magnitude offset to associate events
    self._associate=True
    self._blits = []
    self._drawing=False # are we in a canvas drawing phase?
    self._resizing=False # are we in a window resize phase?
    self._origins = array([]) # dictionary of input origins
    self._triggers = array([]) # dictionary of input origins
    self._refdict = {} # dictionary of reference events
    self._fixmags = True # should we use Sadeh et al., 2014 Pd  - magnitude relations?
    self._fixdist= True # use real distance of estimated distance
    self._deletedEvents = [] # a list of deleted events
    self._background = None
    self._meterBackground = None
    #self.lastEventClicked = None # last element on tree to be clicked
    self.lastEventEntered = None # last element on tree to be entered
    self.auxiliaryLines = [] # a list of auxiliary lines
    self.Bbox = bboxForm()
    self.Bbox.setLims(*self.bbox)
    if splash: splash.showMessage('Initializing...',Qt.AlignCenter)
    QApplication.processEvents()
    QMainWindow.__init__(self, parent)
    self.setWindowTitle('ElarmS Review Tool')
    self.args = args # save command line arguments
    self.meter = mpl.lines.Line2D([],[],color='r') # line for measuring distances along canvas
    if splash: splash.showMessage('Creating application...',Qt.AlignCenter)
    QApplication.processEvents()
    self.create_menu() # create the app top menu
    self.create_status_bar() # add a status bar
    self.create_main_frame() # create the main frame.
    self.create_toolbox() # create a toolbox
    self.create_event_widget() # add events widget.
    self.create_info_widget() # Add info widget.
    self.create_trigger_widget() # Add origin widget.
    self.create_tooltip_widget()
    self.blur = QGraphicsBlurEffect(self.canvas)
    self.canvas.setGraphicsEffect(self.blur)
    self.blur.setEnabled(False)
    self.osm = osm(self.ax,OSMTILEURL,OSMTILEPAT,OSMTILEARCHIVE) # add open street map generator
    self.zoomIsrael()
    if args.i:
      if splash: splash.showMessage('Reading E2 log file(s)...',Qt.AlignCenter)
      QApplication.processEvents()
      self.get_input_file(args.i,0)
    if args.r:
      if splash: splash.showMessage('Reading catalog data...',Qt.AlignCenter)
      QApplication.processEvents()
      self.get_ref_file(args.r,0)
    self.init_connections() # initialize signal connections
    if splash: splash.showMessage('Populating widgets...',Qt.AlignCenter)
    QApplication.processEvents()
    self.refreshDisplay(True)

  def init_connections(self):
    '''Connect signals to functions.
       Using signals to run functions from subprocesses.'''
    self.connect(self,SIGNAL('drawSignal'),self.draw) # draw canvas
    self.connect(self,SIGNAL('blitSignal'),self.blit) # blit draw canvas
    self.Bbox.accepted.connect(self.onBboxAccepted) # connect bbox to dialog ok button
    self.connect(self.minmagLine, SIGNAL('ok'),self.updateMagLims)
    self.connect(self.maxmagLine, SIGNAL('ok'),self.updateMagLims)
    self.starttimeLine.dateTimeChanged.connect(self.updatestarttime)
    self.endtimeLine.dateTimeChanged.connect(self.updateendtime)
    self.fixmagsCheck.stateChanged.connect(self.updateFixMags)
    self.fixdistCheck.stateChanged.connect(self.updateFixDist)
    self.reportedOnlyCheck.stateChanged.connect(self.updateReportedOnly)
    # connecting mpl events
    self.canvas.mpl_connect('scroll_event',self.scroll_event) # connect canvas scroll event
    self.canvas.mpl_connect('after_home_event',self.handle_home) # connect home button on toolbar
    self.canvas.mpl_connect('after_back_event', self.handle_home) # connect back button on toolbar
    self.canvas.mpl_connect('after_forward_event',self.handle_home) # connect forward button on toolbar
    self.canvas.mpl_connect('after_release_pan_event',self.handle_home) # connect pan button on toolbar
    self.canvas.mpl_connect('after_release_zoom_event',self.handle_home) # connect zoom button on toolbar
    self.canvas.mpl_connect('motion_notify_event',self.on_move) # connect mouse motion on canvas
    self.canvas.mpl_connect('button_press_event',self.on_click) # connect click on canvas
    self.canvas.mpl_connect('button_release_event',self.on_unclick) # connect mouse button release on canvas
    self.canvas.mpl_connect('axes_leave_event',self.on_ax_leave) # connect leave map area event on canvas
    #self.canvas.mpl_connect('resize_event',self.resizeEvent) # connect resize event
    # connect eventWidget signals
    self.eventWidget.itemClicked.connect(self.on_event_dbclicked)
    self.eventWidget.itemDoubleClicked.connect(self.on_event_clicked)
    self.eventWidget.itemEntered.connect(self.on_event_entered)
    # connect originWidget signals
    self.historyWidget.itemClicked.connect(self.on_history_dbclicked)
    self.historyWidget.itemDoubleClicked.connect(self.on_history_clicked)
    self.historyWidget.itemEntered.connect(self.on_history_entered)
    # connect magnitude histogram and plot events
    self.magplotWidget.ax.figure.canvas.mpl_connect('pick_event',self.on_line_Pick) # pick event

  def draw(self,idle=True):
    'Draw the figure on canvas.'
    if self._drawing: return
    self._drawing = True # entering a drawing mode
    if self.lastEventEntered: self.event_restore(self.lastEventEntered)
    if idle:
      self.canvas.draw_idle() # draw idle (see matplotlib for details)
    else:
      self.canvas.draw() # draw anyways.
    self._drawing = False # leaving drawing mode
    QApplication.processEvents()
    self._background = self.canvas.copy_from_bbox(self.ax.bbox)

  def blit(self):
    self.canvas.restore_region(self._background)
    [self.ax.draw_artist(a) for a in self._blits]
    self._blits=[]
    self.canvas.blit(self.ax.bbox)

  def showEvent(self,event):
    self.resizeEvent(None)
    self._background = self.canvas.copy_from_bbox(self.ax.bbox)

  def resizeEvent(self,event):
    '''called by any resize event of the map'''
    if not self._resizing: # only run if we are not in a middle of a resizong process
      self._resizing=True # note that we are resizong
      self.redrawbgmap() # redraw the background map
      self._resizing=False # Done with resizing.
#    event.accept()

  def zoom(self,e):
    '''called by scroll_event to zoom in or out on map.'''
    ax = self.ax # easier
    x,y = e.xdata,e.ydata # get where to center to (mouse pointer)
    x0,x1 = ax.get_xlim() # get current x limits
    y0,y1 = ax.get_ylim() # get current y limits
    dx=(x1-x0)/2.0 # get distance from center to edge on x axis
    dy=(y1-y0)/2.0 # get distance from center to edge on y axis
    # zoom and center
    ax.set_xlim(x-dx*e.zoom,x+dx*e.zoom)
    ax.set_ylim(y-dy*e.zoom,y+dy*e.zoom)
    self.redrawbgmap() # redraw background map with new limits

  def scroll_event(self,e):
    '''called by a scroll event on map'''
    if e.button=='up': e.zoom = 0.9 # zoom in
    elif e.button=='down': e.zoom = 1.1 # zoom up
    self.zoom(e) # do the zoom

  def redrawbgmap(self):
    '''redraws the background map using open street map object (osm)
    see osm module for more details.
    '''
    self.statusBar().showMessage('hold, redrawing ...') # let user know to wait on drawing. it might take some time
    self.blur.setEnabled(True)
    QApplication.processEvents() # make sure user see the message
    self.ax.images=[] # remove images from map
    self.ax.apply_aspect() # make sure xlim and ylim are updated to screen size. this is because we use equal aspect and datalim. see matplotlib details on axes set_aspect function
    x0,x1 = self.ax.get_xlim() # get requested limits of x axis
    y0,y1 = self.ax.get_ylim() # get requested limits of y axis
    self.osm.relimcorrected(x0,x1,y0,y1) # make sure limits are not out of map phisical boundaries. see osm module for more details.
    tiles = self.osm.gettiles(x0,x1,y0,y1) # get the needed tiles from buffer or url. see osm module for more details
    self.osm.plottiles(tiles) # plot the tiles. see osm module for more details
    self.emit(SIGNAL('drawSignal')) # call self draw using a signal in case we are in a subprocess.
    self.statusBar().showMessage('Done redrawing.',2000) # note user we are done with redrawing.
    self.blur.setEnabled(False)
    self.blur.update()
    QApplication.processEvents() # make sure user see the message

  def handle_home(self,evt):
    'handle mpl toolbar pan/zoom/back/forward/home buttons.'
    self.fixlimits() # fix axes limit to global geographic bounds
    self.redrawbgmap() # redraw the background map

  def on_move(self,evt):
    'handle mouse movements along map'
    if evt.inaxes and self.toolbar.mode=='': # make sure we are not in a toolbar mode of pan/zoom
      if evt.button==3: # check if we measure distance along the map (right button is clicked)
        self.meter.set_data([self.meter._xorig[0],evt.xdata],[self.meter._yorig[0],evt.ydata]) # adjust meter line edges
        self._blits.append(self.meter)
        self.statusBar().showMessage('Distance: %lfkm'%(geo.inv(self.meter._xorig[0],self.meter._yorig[0],self.meter._xorig[-1],self.meter._yorig[-1])[-1]/1000.)) # show user the distance
        self.emit(SIGNAL('blitSignal'),True) # draw the map
        return # we're done here
      hide = True # unless we don't measure but simply pointing along the map
      label = [] # a list of labels to present in a tooltip
      for l in self.ax.hitlist(evt): # see if mouse points at objects of the map
        if l in self.ax.lines and not l in self.auxiliaryLines: # if object is a line
          if not l.get_label() in label: label.append(l.get_label()) # get it's label to the list
          hide=False # a flag for tooltip widget
      if not hide: # if we don't hide the tooltip
        self.stationNameWidget.setText('\n---\n'.join(label)) # set text in tooltip
        self.stationNameWidget.move(QCursor.pos()+QPoint(5,5)) # move the tooltip to the pointer position. might be out of screen area if on edges
        self.stationNameWidget.adjustSize() # adjust the size of tooltip widget to the text
        self.stationNameWidget.show() # show the widget
      else: self.stationNameWidget.hide() # if nothing is on the hit list - hide the tooltip widget

  def on_click(self,evt):
    'handle mouse clicks'
    if evt.button==1 and not evt.dblclick: # in case its a single left button click
      if evt.inaxes and self.toolbar.mode=='': # and not on a toolbar action of pan/zoom
        for l in self.ax.hitlist(evt): # see if mouse points at objects of the map
          if type(l)==mpl.lines.Line2D:
            try:
              item = l.item # get the item
            except:
              continue
            if type(item)==StationMagnitude: continue
            if type(item)==Origin: item=item.event
            if not type(item)==CatalogEvent: continue
            self.eventWidget.clearSelection()# unselect other items
            self.eventWidget.setItemSelected(item,True) # select the item
            self.eventWidget.scrollToItem(item) # go to item on list
            self.on_event_dbclicked(item)
            return # we're out of here
    elif evt.button==3 and evt.dblclick: # if its a right button double click
        return # we're out of here
    elif evt.button==3 and not evt.dblclick: # if its a right button and not a double click
      if evt.inaxes and self.toolbar.mode=='': # and not on a toolbar action of pan/zoom
        self.meter.set_data([evt.xdata],[evt.ydata]) # update the meter measurement tool
        self._meterBackground = self._background
        self.ax.add_artist(self.meter) # add the meter to the axes
        self._blits.append(self.meter) # add the meter to the blits list
        self.emit(SIGNAL('blitSignal')) # fast drawing
        return # we're out of here

  def on_unclick(self,evt):
    'handle mouse unclick or mouse button release'
    if evt.button==3 and any(self.meter.get_data()): # if its the right button and there is some data set in the meter
      self.meter.remove() # remove the meter from the axes
      self._background = self._meterBackground
      self._meterBackground = None
      self.emit(SIGNAL('blitSignal'),True) # redraw the map
      self.statusBar().showMessage('Distance: %lfkm'%(geo.inv(self.meter._xorig[0],self.meter._yorig[0],self.meter._xorig[-1],self.meter._yorig[-1])[-1]/1000.),1000) # last update of the measurement distance
      self.meter.set_data([],[]) # remove any data from the meter line

  def on_ax_leave(self,event):
    'things to do when leaving the map area with mouse'
    self.stationNameWidget.hide() # hide the tooltip widget

  def keyPressEvent(self,event):
    'Handle key pressed'
    if event.isAutoRepeat(): return
    if event.key()==Qt.Key_Escape:
      self.clearSelectedEvents()
    if event.key()==Qt.Key_Delete:
      items = self.eventWidget.selectedItems()
      if not len(items): return
      for item in items:
        self.deletEvents(item)
    if event.key()==Qt.Key_Z and (event.modifiers() & Qt.ControlModifier):
      self.undeleteEvents()

  def clearSelectedEvents(self):
      self.auxiliaryClear()
      if self.lastEventEntered:
        self.event_restore(self.lastEventEntered)
        self.lastEventEntered = None
      self.eventWidget.clearSelection()
      self.historyWidget.clearSelection()
      self.triggerWidget.clearSelection()
      self.emit(SIGNAL("drawSignal"))

  def grid(self):
    'toggle grid on or off'
    self.ax.grid()
    self.draw()

  def create_tooltip_widget(self):
    'creates tooltip of figure elements'
    self.stationNameWidget = QLabel() # take a QLabel
    self.stationNameWidget.setFrameShape(QFrame.StyledPanel) # add a frame
    self.stationNameWidget.setWindowFlags(Qt.ToolTip) # make window look like a tooltip
    self.stationNameWidget.setAttribute(Qt.WA_TransparentForMouseEvents) # mouse events can't affect it.
    self.stationNameWidget.hide() # hide for now. see self.on_move function on how to use.

  def drag_pan(self,button,key,x,y):
    'a replacement to the original drag_pan function of the mpl axes.'
    mpl.axes.Axes.drag_pan(self.ax,button,key,x,y) # see matplotlib for more details.
    self.fixlimits() # fix the limits so we stay within the geo bounderies.

  def fixlimits(self):
    'fix the limits of the axes to acceptable geographic bounds'
    x0,x1 = self.ax.get_xlim() # get current longitude limits
    y0,y1 = self.ax.get_ylim() # get current latitude limits
    dx = x1-x0 # distance from center to edge horizontal
    dy = y1-y0 # distance from center to edge vertical
    if dy>150: dy=150 # don't exceed latitude limit
    if y1>75: # fix maximal latitude
      y1=75
      y0=y1-dy # and adjust minimal one
    if y0<-75: # fix maximal latitude
      y0 = -75
      y1 = y0+dy # and adjust maximal one
    if dx>358: dx=358 # don't exceed longidute limit
    if x1>179: # fix maximal longitude
      x1=179
      x0=x1-dx # adjust minimal one
    if x0<-179: # fix minimal longitude
      x0=-179
      x1=x0+dx# adjust maximal one
    # set limits to corrected ones
    self.ax.set_ylim(y0,y1)
    self.ax.set_xlim(x0,x1)

  def goToLocation(self,l):
    'zoom to a location'
    x,y = (l._x[0],l._y[0]) # get the location
    self.osm.relimcorrected(x-0.5,x+0.5,y-0.5,y+0.5) # change limits of the map to zoom in on home location
    self.redrawbgmap() # redraw the map (updating background maps)

  def zoomIsrael(self):
    'zoom to israel'
    self.osm.relimcorrected(32.5,37.5,29,34) # adjust map limits
    self.redrawbgmap() # redraw the map (updating background maps)

  def create_event_widget(self):
    'create the event list widget'
    self.eventWidget = QTreeWidget()
    self.eventWidget.setAlternatingRowColors(True)
    self.eventWidget.setColumnCount(13)
    self.eventWidget.setSortingEnabled(True)
    self.eventWidget.setMouseTracking(True)
    header = ['EvID','Origin-Time','Mag','Latitude','Longitude','Origin-Time','Mag','Latitude','Longitude','Loc Err','OT Err','RefID','Event Type']
    self.eventWidget.setHeaderLabels(header)
    self.eventWidget.sortItems(1,1)
    self.treeSplitter.addWidget(self.eventWidget)
    self.treeSplitter.setStretchFactor(1,1)

  def create_history_widget(self):
    'creates the event history widget.'
    self.historyWidget = QTreeWidget()
    self.historyWidget.setAlternatingRowColors(True)
    self.historyWidget.setColumnCount(13)
    header = ['EvID','Alert-Time','Origin-Time','Mag','Latitude','Longitude','Loc Err','OT Err','nT','nS','%S','Alert','OK']
    self.historyWidget.setHeaderLabels(header)
    self.historyWidget.setSortingEnabled(True)
    self.historyWidget.setMouseTracking(True)
    self.historyWidget.sortItems(1,0)

  def create_trigger_widget(self):
    'create the triggers list widget'
    self.triggerWidget = QTreeWidget()
    self.triggerWidget.setAlternatingRowColors(True)
    self.triggerWidget.setColumnCount(16)
    self.triggerWidget.setSortingEnabled(True)
    self.triggerWidget.setMouseTracking(True)
    header = ['EvID','StaID','Trigger-Time','Latitude','Longitude','DistKm','TTErr','Mag (Pd)','Used','LogPd','PdSNR','Mag (tauP)','Used','LogTauP','TaupSNR','Azimuth']
    self.triggerWidget.setHeaderLabels(header)
    self.triggerWidget.sortItems(2,0)
    self.treeSplitter.addWidget(self.triggerWidget)
    self.treeSplitter.setStretchFactor(3,1)

  def create_magerr_widget(self):
    self.magerrWidget = QWidget()
    fig = Figure(figsize=(1,1),dpi=100)
    canvas = FigureCanvas(fig)
    self.magerrWidget.ax = fig.add_subplot(111)
    self.magerrWidget.ax.plot([],[],'kv',mfc='None')
    self.magerrWidget.ax.set_ylabel('Magnitude Error (E2 - Catalog)')
    self.magerrWidget.ax.set_xlabel('Time since origin time (Alert - Catalog)')
    v = QVBoxLayout(self.magerrWidget)
    tb = NavigationToolbar(canvas,self.magerrWidget)
    canvas.setParent(self.magerrWidget)
    v.addWidget(canvas)
    v.addWidget(tb)

  def create_timeerr_widget(self):
    self.timeerrWidget = QWidget()
    fig = Figure(figsize=(1,1),dpi=100)
    canvas = FigureCanvas(fig)
    self.timeerrWidget.ax = fig.add_subplot(111)
    self.timeerrWidget.ax.plot([],[],'kv',mfc='None')
    self.timeerrWidget.ax.set_ylabel('Origin Time Error (E2 - Catalog)')
    self.timeerrWidget.ax.set_xlabel('Time since origin time (Alert - Catalog)')
    v = QVBoxLayout(self.timeerrWidget)
    tb = NavigationToolbar(canvas,self.timeerrWidget)
    canvas.setParent(self.timeerrWidget)
    v.addWidget(canvas)
    v.addWidget(tb)

  def create_locerr_widget(self):
    self.locerrWidget = QWidget()
    fig = Figure(figsize=(1,1),dpi=100)
    canvas = FigureCanvas(fig)
    self.locerrWidget.ax = fig.add_subplot(111)
    self.locerrWidget.ax.plot([],[],'kv',mfc='None')
    self.locerrWidget.ax.set_ylabel('Location Error (km)')
    self.locerrWidget.ax.set_xlabel('Time since origin time (Alert - Catalog)')
    v = QVBoxLayout(self.locerrWidget)
    tb = NavigationToolbar(canvas,self.locerrWidget)
    canvas.setParent(self.locerrWidget)
    v.addWidget(canvas)
    v.addWidget(tb)

  def create_summary_widget(self):
    self.summaryWidget = QTableWidget(5,2)
    self.summaryWidget.verticalHeader().setHidden(True)
    self.summaryWidget.horizontalHeader().setHidden(True)
    self.summaryWidget.setGridStyle(Qt.NoPen)
    self.summaryWidget.setAlternatingRowColors(True)
    self.summaryWidget.setSelectionMode(QAbstractItemView.NoSelection)
    self.summaryWidget.setItem(0,0,QTableWidgetItem("Events"))
    self.summaryWidget.setItem(1,0,QTableWidgetItem("Matched"))
    self.summaryWidget.setItem(2,0,QTableWidgetItem("Missed"))
    self.summaryWidget.setItem(3,0,QTableWidgetItem("False"))
    self.summaryWidget.setItem(4,0,QTableWidgetItem("Score"))
    self.summaryWidget.setItem(0,1,QTableWidgetItem("    -"))
    self.summaryWidget.setItem(1,1,QTableWidgetItem("    -"))
    self.summaryWidget.setItem(2,1,QTableWidgetItem("    -"))
    self.summaryWidget.setItem(3,1,QTableWidgetItem("    -"))
    self.summaryWidget.setItem(4,1,QTableWidgetItem("    -"))

  def create_maghist_widget(self):
    self.maghistWidget =  QWidget()
    fig = Figure(figsize=(1,1),dpi=100)
    canvas = FigureCanvas(fig)
    self.maghistWidget.ax = fig.add_subplot(111)
    self.maghistWidget.ax.set_ylabel('# Magnitude Solutions')
    self.maghistWidget.ax.set_xlabel('Magnitude Error (ElarmS - Catalog)')
    v = QVBoxLayout(self.maghistWidget)
    tb = NavigationToolbar(canvas,self.maghistWidget)
    canvas.setParent(self.maghistWidget)
    v.addWidget(canvas)
    v.addWidget(tb)

  def create_magplot_widget(self):
    self.magplotWidget = QWidget()
    fig = Figure(figsize=(1,1),dpi=100)
    canvas = FigureCanvas(fig)
    self.magplotWidget.ax = fig.add_subplot(111)
    self.magplotWidget.ax.set_ylabel('ElarmS Magnitude')
    self.magplotWidget.ax.set_xlabel('Catalog Magnitude')
    v = QVBoxLayout(self.magplotWidget)
    tb = NavigationToolbar(canvas,self.magplotWidget)
    canvas.setParent(self.magplotWidget)
    v.addWidget(canvas)
    v.addWidget(tb)

  def create_info_widget(self):
    'creates the info widget.'
    self.infoWidget = QTabWidget(self.treeSplitter)
    self.treeSplitter.setStretchFactor(2,1)
    self.create_summary_widget()
    self.infoWidget.addTab(self.summaryWidget,'Summary')
    self.create_maghist_widget()
    self.infoWidget.addTab(self.maghistWidget,'Histogram')
    self.create_magplot_widget()
    self.infoWidget.addTab(self.magplotWidget,'Plot')
    self.create_history_widget()
    self.infoWidget.addTab(self.historyWidget,'History')
    self.create_timeerr_widget()
    self.infoWidget.addTab(self.timeerrWidget,'Time Error')
    self.create_magerr_widget()
    self.infoWidget.addTab(self.magerrWidget,'Mag Error')
    self.create_locerr_widget()
    self.infoWidget.addTab(self.locerrWidget,'Loc Error')

  def save_preffered_triggeres(self):
    'Save parameters of event and triggers'
    data = [['cEID','cOT','cMAG','cLON','cLAT',
             'oEID','oOT','oMAG','oLON','oLAT',
             'tID','tT','tLON','tLAT','logpd','pdmag','logtp','dist']]
    for i in xrange(self.eventWidget.topLevelItemCount()):
      ev = self.eventWidget.topLevelItem(i)
      if ev.isHidden(): continue # skip hidden
      if ev.preffered: 
        origin = ev.preffered # get preffered origin
      else:
        continue # don't waste time on no preffered events
      for trig in origin.triggers:
        data += [[ev.EID,str(ev.ot),str(ev.mag),str(ev.lon),str(ev.lat),
             origin.EID,str(origin.ot),str(origin.mag),str(origin.lon),str(origin.lat),
             trig.ID,str(trig.tt),str(trig.lon),str(trig.lat),str(trig.logpd),str(trig.pdmag),str(trig.logtp),str(trig.dist)]]
    fileurl = str(QFileDialog.getSaveFileName(self, 'Save Table Data',
                                                filter='CSV Files (*.csv);;Space Separated Files (*.txt);;All Files (*.*)')) # get the file name
    if not fileurl: return
    if fileurl.endswith('txt'):
      sep = ' '
    else:
      sep = ','
    try:
      with open(fileurl,'w') as f:
        for d in data:
          f.write(sep.join(d)+'\n')
    except:
      self.message('Saving data to %s failed'%fileurl)
      return
    self.statusBar().showMessage('Saved %d events to %s'%(len(data),fileurl), 2000)

  def get_input_file(self,filesurl=None,update=1):
    'open a dialog to get a file name'
    if not filesurl:
      filesurl = QFileDialog.getOpenFileNames(self, 'Open ElarmS Log file(s)',filter='Log Files (*.log);;All Files (*.*)') # get the file name
    LOGs = [str(f) for f in filesurl]
    self._origins,self._triggers = concatenateLOG(LOGs)
    if not len(self._origins):
      self.statusBar().showMessage('No Events in file.', 2000)
      if update: self.refreshDisplay(True)
      return
    self.statusBar().showMessage('Loaded %d events'%len(self._origins), 2000)
    self._associate = True
    if update: self.refreshDisplay(True)

  def get_ref_file(self,filesurl=None,update=1):
    'open a dialog to get a file name'
    if not filesurl:
      filesurl = QFileDialog.getOpenFileNames(self, 'Open Reference data file(s)',filter='CSV Files (*.csv);;All Files (*.*)') # get the file name
    CSVs = [str(f) for f in filesurl]
    self._refdict = concatenateCSV(CSVs)
    if not len(self._refdict):
      self.statusBar().showMessage('No Events in file.', 2000)
      if update: self.refreshDisplay(True)
      return
    self.statusBar().showMessage('Loaded %d events'%len(self._refdict), 2000)
    self._associate = True
    if update: self.refreshDisplay(True)

  def updatetables(self):
    self.clearTables()
    print UTCDateTime(),'updatetables'
    if self._associate:
      print UTCDateTime(),'updatetables associate'
      self.associate()
    # get associated and unassociated events
    refs = [v['refID'] for v in self._refdict.values()] #  associated
    unassociated = dict([(EID,'') for EID in self._origins if not EID in refs]) # unassociated
    # scan catalog events
    print UTCDateTime(),'updatetables associated'
    for ref in self._refdict.values():
      parent = CatalogEvent(*ref) # create an event line
      if parent.refID in self._origins: # look for associated origins
        parent.l.set_color('yellow') # mark as True
        for origin in self._origins[parent.refID]:# look for all origins
          event = Origin(*origin) # create an origin line
          parent.addOrigin(event) # add origin to event
          event.l.set_color('green') # mark as True
          if event.first and not event.l in self.ax.lines: self.ax.add_line(event.l)# add line to map
          self.historyWidget.addTopLevelItem(event)# add origin to history widget
          event.setHidden(True)
          for trigger in self._triggers['-'.join([event.EID,str(event.ver)])]:
            trig = StationMagnitude(*trigger)
            if not trig.l in self.ax.lines: self.ax.add_line(trig.l)# add line to map
            event.addTrigger(trig)
            self.triggerWidget.addTopLevelItem(trig)
            trig.setHidden(True)
      else:
        parent.l.set_color('orange') # mark as Missed
      if not parent.l in self.ax.lines: self.ax.add_line(parent.l)# add line to map
      self.eventWidget.addTopLevelItem(parent) # add item to widget
    print UTCDateTime(),'updatetables unassociated'
    # scan unassociated origins
    for EID in unassociated:
      parent = CatalogEvent(refID=EID) # create an event line
      if not len([origin for origin in self._origins[EID] if origin['first']])>0: continue # don't add unassociated unreported events
      for origin in self._origins[EID]: # look for all origins
        event = Origin(*origin) # create an origin line
        parent.addOrigin(event) # add origin to event
        event.l.set_color('red') # mark as False
        if event.first:
          if not event.l in self.ax.lines: self.ax.add_line(event.l)# add line to map
          parent.ot = event.ot # add origin time to event
          parent.setText(1,str(event.ot))
        self.historyWidget.addTopLevelItem(event)# add origin to history widget
        event.setHidden(True)
        for trigger in self._triggers['-'.join([event.EID,str(event.ver)])]:
          trig = StationMagnitude(*trigger)
          if not trig.l in self.ax.lines: self.ax.add_line(trig.l)# add line to map
          event.addTrigger(trig)
          self.triggerWidget.addTopLevelItem(trig)
          trig.setHidden(True)
      if not parent.l in self.ax.lines: self.ax.add_line(parent.l)# add line to map
      self.eventWidget.addTopLevelItem(parent) # add item to widget

  def refreshDisplay(self,update=False):
    self.blur.setEnabled(True)
    QApplication.processEvents()
    if update:
      self.updatetables()
    print UTCDateTime(),'updatetables Fix Mag'
    self.fixMags()
    print UTCDateTime(),'updatetables Filter'
    self.filterEvents() # hide unwanted events
    print UTCDateTime(),'updatetables summary'
    self.updateSummary()
    print UTCDateTime(),'updatetables histogram and plot'
    self.updateHistPlot()
    print UTCDateTime(),'updatetables done'
    self.emit(SIGNAL('drawSignal'))
    self.refreshbuttonHilight.setStrength(0) # reset button colors
    self.blur.setEnabled(False)
    self.blur.update()
    QApplication.processEvents()

  def clearTables(self):
    # remove points from plot
    self.clearMap()
    # clear tables
    self.eventWidget.clear()
    self.historyWidget.clear()
    self.triggerWidget.clear()
    # clear plots
    self.timeerrWidget.ax.lines[0].set_data([],[])
    self.timeerrWidget.ax.figure.canvas.draw()
    self.timeerrWidget.ax.lines[0].set_data([],[])
    self.timeerrWidget.ax.figure.canvas.draw()
    self.timeerrWidget.ax.lines[0].set_data([],[])
    self.timeerrWidget.ax.figure.canvas.draw()

  def clearMap(self):
    # remove points from plot
    for i in xrange(self.eventWidget.topLevelItemCount()):
      event = self.eventWidget.topLevelItem(i)
      for origin in event.origins:
        for trigger in origin.triggers:
          if trigger.l in self.ax.lines: self.ax.lines.remove(trigger.l)
        if origin.l in self.ax.lines: self.ax.lines.remove(origin.l)
      if event.l in self.ax.lines: self.ax.lines.remove(event.l)


  def filterEvents(self):
    'filter events by time,location and magnitude limits'
    for i in xrange(self.eventWidget.topLevelItemCount()):
      FILTER=False # don't filter default
      ev = self.eventWidget.topLevelItem(i) # get an event from the list
      event = ev # generate a reference
      if not ev.lat or not ev.lon: # if its an unassociated event
        event = ev.preffered # get its prefered origin as reference
      if not (self.minmag <= event.mag <= self.maxmag \
              and UTCDateTime(self.starttime) <= event.ot <= UTCDateTime(self.endtime) \
              and self.bbox[0] <= event.lon <=self.bbox[1] \
              and self.bbox[2] <= event.lat <= self.bbox[3]):
        FILTER=True # filter if the reference is out of limits
      try:
        if FILTER and (self.minmag <= event.preffered.mag <= self.maxmag \
              and UTCDateTime(self.starttime) <= event.preffered.ot <= UTCDateTime(self.endtime) \
              and self.bbox[0] <= event.preffered.lon <=self.bbox[1] \
              and self.bbox[2] <= event.preffered.lat <= self.bbox[3]):
          FILTER=False # don't filter if the reference preffered origin is within limits
      except: # we are already looking at an un associated origin
        pass
      if FILTER:
        ev.setHidden(True)
        if ev.l in self.ax.lines: self.ax.lines.remove(ev.l)
        try:
          if ev.preffered.l in self.ax.lines: self.ax.lines.remove(ev.preffered.l)
        except:
          pass
      else:
        ev.setHidden(False)
        if not ev.l in self.ax.lines: self.ax.add_line(ev.l)
        try:
          if not ev.preffered.l in self.ax.lines: self.ax.add_line(ev.preffered.l)
        except:
          pass

  def filterBlasts(self):
    for i in xrange(self.eventWidget.topLevelItemCount()):
      event = self.eventWidget.topLevelItem(i)
      if int(event.text(12))==2:
        event.setHidden(True)
        if event.l in self.ax.lines: self.ax.lines.remove(event.l)
        try:
          if event.preffered.l in self.ax.lines: self.ax.lines.remove(event.preffered.l)
        except:
          pass
    self.updateSummary()
    self.updateHistPlot()
    self.draw()

  def deletEvents(self,item):
    item.setHidden(True)
    if item.l in self.ax.lines: self.ax.lines.remove(item.l)
    try:
      if item.preffered.l in self.ax.lines: self.ax.lines.remove(item.preffered.l)
    except:
      pass
    self.updateSummary()
    self.updateHistPlot()
    self.clearSelectedEvents()
    self._deletedEvents.append(item)
    self.draw()

  def undeleteEvents(self):
    if not len(self._deletedEvents): return
    item = self._deletedEvents.pop(-1)
    item.setHidden(False)
    if not item.l in self.ax.lines: self.ax.add_line(item.l)
    try:
      if not item.preffered.l in self.ax.lines: self.ax.add_line(item.preffered.l)
    except:
      pass
    self.updateSummary()
    self.updateHistPlot()
    self.clearSelectedEvents()
    self.draw()

  def associate(self):
    self._associate=False
    assoc = []
    for event in self._refdict.values():
      bestOrigin = None
      for refID in self._origins:
        if refID in assoc: continue
        origin = [origin for origin in self._origins[refID] if origin['first']] # get the first reported origin
        if not len(origin)>0: # if non was reported...
          if self.reportedOnly: continue # if we use only reported events skip this one
          origin = [orig for orig in self._origins[refID] if orig['OK']==7] # get all origins that didn't alert due to low magnitude
          if not len(origin)>0: # if no origin is available...
            continue
        origin = origin[0] # convert list to object
        origin['first']=1 # set as first
        # chack distance and time
        if geo.inv(origin['lon'],origin['lat'],event['lon'],event['lat'])[-1]/1000.>self.deltaR \
           or abs(UTCDateTime(origin['ot'])-UTCDateTime(event['ot']))>self.deltaT: continue
        if self.deltaM and abs(origin['mag']-event['mag'])>self.deltaM: continue
        if not bestOrigin or \
               abs(UTCDateTime(bestOrigin['ot'])-UTCDateTime(event['ot']))>abs(UTCDateTime(origin['ot'])-UTCDateTime(event['ot'])) or \
               geo.inv(bestOrigin['lon'],origin['lat'],event['lon'],event['lat'])[-1]/1000.>geo.inv(origin['lon'],origin['lat'],event['lon'],event['lat'])[-1]/1000.:
          bestOrigin=origin
      if bestOrigin:
        event['refID']=bestOrigin['EID']
        assoc.append(bestOrigin['EID'])

  def fixMags(self):
    'using Sadeh et al., 2014 to fix magnitudes'
    for i in xrange(self.eventWidget.topLevelItemCount()):
      event = self.eventWidget.topLevelItem(i)
      if not len(event.origins): continue
      for origin in event.origins:
        mags = []
        for trigger in origin.triggers:
          if self._fixdist:
            dist=geo.inv(trigger.lon,trigger.lat,event.lon,event.lat)[-1]/1000.
          else:
            dist=trigger.dist
          if self._fixmags:
            trigger.pdmag = Sadehpdmag(trigger.logpd,dist)
          else:
            trigger.pdmag = E2pdmag(trigger.logpd,dist)
          if trigger.updm: mags.append(trigger.pdmag)
          trigger.setText(7,"%.2f"%trigger.pdmag)
        if not len(mags): continue
        origin.mag = array(mags).mean()
        origin.setText(3,"%.2f"%origin.mag)
      if event.preffered:
        event.setText(6,"%.2f"%event.preffered.mag)

  def on_event_entered(self,item,col=0):
    'event tree element is entered by mouse'
    if self.lastEventEntered==item: return
    if self.lastEventEntered: self.event_restore(self.lastEventEntered)
    self.emit(SIGNAL("blitSignal"))
    self._background = self.canvas.copy_from_bbox(self.ax.bbox)
    item.l.set_markersize(10) # enlarge
    item.l.set_zorder(99) # bring to front
    self._blits.append(item.l)
    if item.preffered:
      item.preffered.l.set_markersize(10) # enlarge
      item.preffered.l.set_zorder(99) # bring to front
      self._blits.append(item.preffered.l)
    self.lastEventEntered = item
    self.emit(SIGNAL("blitSignal"))

  def on_history_entered(self,item,col=0):
    self.on_event_entered(item.event,0)

  def on_event_clicked(self,item,col=0):
    'click on event widget element'
    if item.lat and item.lon:
      self.goToLocation(item.l) # zoom to item
    else:
      self.goToLocation(item.preffered.l)
    #self.emit(SIGNAL("drawSignal"))

  def on_history_clicked(self,item,col=0):
    'click on event widget element'
    self.goToLocation(item.l) # zoom to item

  def on_event_dbclicked(self,item,col=0):
    self.updateHistory(item)
    self.updateInfo(item)
    if item.preffered:
      self.historyWidget.setItemSelected(item.preffered,True)
      self.historyWidget.scrollToItem(item.preffered)
      self.on_history_dbclicked(item.preffered)
    #if not self._drawing: self.emit(SIGNAL("drawSignal"))

  def on_history_dbclicked(self,item,col=0):
    self.updateTriggers(item)
    self.auxiliarySet(item.event,item)
    if not item==item.event.preffered: # show the event location if not the preffered one
      if not item.l in self.ax.lines:
        self.ax.add_line(item.l)
        self.auxiliaryLines.append(item.l)
    if not self._drawing: self.emit(SIGNAL("drawSignal"))

  def on_line_Pick(self,event):
    try:
      event.canvas.toolbar.set_message(event.artist.get_label())
    except Exception as E:
      print str(E)

  def updateInfo(self,item):
    if item.lat and item.lon:
      x = [] # time since origin
      t = [] # ot err
      m = [] # mag err
      l = [] # loc err
      reported = [] # reported origin
      for origin in item.origins:
        x.append(origin.AT-item.ot)
        t.append(origin.Secs)
        m.append(origin.mag-item.mag)
        l.append(origin.km)
        reported.append(origin.Aok)
      x,t,m,l = array(x),array(t),array(m),array(l)
      self.timeerrWidget.ax.lines[0].set_data(x,t)
      self.timeerrWidget.ax.collections = [self.timeerrWidget.ax.scatter(x[where(reported)[0]],t[where(reported)[0]],edgecolor='b',marker='o',s=120,facecolor='None')]
      self.magerrWidget.ax.lines[0].set_data(x,m)
      self.magerrWidget.ax.collections = [self.magerrWidget.ax.scatter(x[where(reported)[0]],m[where(reported)[0]],edgecolor='b',marker='o',s=120,facecolor='None')]
      self.locerrWidget.ax.lines[0].set_data(x,l)
      self.locerrWidget.ax.collections = [self.locerrWidget.ax.scatter(x[where(reported)[0]],l[where(reported)[0]],edgecolor='b',marker='o',s=120,facecolor='None')]
    [ax.relim() for ax in [self.timeerrWidget.ax,self.magerrWidget.ax,self.locerrWidget.ax]]
    [ax.autoscale() for ax in [self.timeerrWidget.ax,self.magerrWidget.ax,self.locerrWidget.ax]]
    [ax.figure.canvas.draw() for ax in [self.timeerrWidget.ax,self.magerrWidget.ax,self.locerrWidget.ax]]

  def auxiliarySet(self,item,preffered):
    self.auxiliaryClear()
    if item.lat and item.lon:
      self.auxiliaryLines.append(mpl.lines.Line2D([preffered.lon,item.lon],[preffered.lat,item.lat],color='0.5'))
    for trigger in preffered.triggers:
      self.auxiliaryLines.append(mpl.lines.Line2D([preffered.lon,trigger.lon],[preffered.lat,trigger.lat],linestyle='dotted',color='0.5'))
    [self.ax.add_line(l) for l in self.auxiliaryLines]

  def updateHistory(self,item):
    self.historyWidgetHideAll()
    for origin in item.origins:
      origin.setHidden(False)

  def updateTriggers(self,preffered):
    self.triggerWidgetHideAll()
    for trigger in preffered.triggers:
      trigger.setHidden(False)

  def historyWidgetHideAll(self):
    self.triggerWidgetHideAll()
    if not self.historyWidget.itemAt(0,0): return
    event = self.historyWidget.itemAt(0,0).event
    for origin in event.origins:
      origin.setHidden(True)

  def triggerWidgetHideAll(self):
    if not self.triggerWidget.itemAt(0,0): return
    origin = self.triggerWidget.itemAt(0,0).origin
    for trigger in origin.triggers:
      trigger.setHidden(True)

  def event_restore(self,item):
    'unclick tree element'
    item.l.set_markersize(6) # resize
    item.l.set_zorder(1) # reorder
    self._blits.append(item.l)
    if item.preffered:
      item.preffered.l.set_markersize(6) # enlarge
      item.preffered.l.set_zorder(1) # bring to front
      self._blits.append(item.preffered.l)

  def auxiliaryClear(self):
    [self.ax.lines.remove(l) for l in self.auxiliaryLines]
    self.auxiliaryLines = []

  def updateSummary(self):
    tb = self.summaryWidget
    t,m,f = 0,0,0
    item = None
    # got to first visible item
    for i in xrange(self.eventWidget.topLevelItemCount()):
      item = self.eventWidget.topLevelItem(i)
      if not item.isHidden(): break
    while item:
      if item.lat and item.lon:
        if len(item.origins):
          t=t+1
        else:
          m=m+1
      else:
        f=f+1
      item= self.eventWidget.itemBelow(item)
    if t+m+f:
      tb.item(0,1).setText(str(t+m+f))
      tb.item(1,1).setText(str(t))
      tb.item(2,1).setText(str(m))
      tb.item(3,1).setText(str(f))
      tb.item(4,1).setText("%.2f%%"%(100.0*t/float(t+m+f)))
    else:
      tb.item(0,1).setText('    -')
      tb.item(1,1).setText('    -')
      tb.item(2,1).setText('    -')
      tb.item(3,1).setText('    -')
      tb.item(4,1).setText('    -')


  def updateHistPlot(self):
    if not self.eventWidget.topLevelItemCount(): return
    tmag = []
    tmagu = []
    omag = []
    for i in xrange(self.eventWidget.topLevelItemCount()):
      event = self.eventWidget.topLevelItem(i)
      if not event.lat and not event.lon: continue
      if event.isHidden(): continue
      if not event.preffered: continue
      omag.append((event.mag,event.preffered.mag,event.refID))
      for trigger in event.preffered.triggers:
        if trigger.updm:
          tmag.append((event.mag,trigger.pdmag,event.refID+' '+trigger.ID))
        else:
          tmagu.append((event.mag,trigger.pdmag,event.refID+' '+trigger.ID))
    ax = self.maghistWidget.ax
    ax.patches=[]
    ax.lines=[]
    if len(tmag) and len(omag):
      dm = array([(mag[1]-mag[0]) for mag in tmag])
      ax.hist(dm,50,(-5,5),histtype='stepfilled',color='0.5',label='Stations (%d) mean: %.2f std: %.2f'%(len(dm),dm.mean(),dm.std()))
      dm = array([(mag[1]-mag[0]) for mag in omag])
      ax.hist(dm,50,(-5,5),histtype='stepfilled',color='0.25',label='Events (%d) mean: %.2f std: %.2f'%(len(dm),dm.mean(),dm.std()))
      ax.legend(frameon=False,loc=6,handlelength=0.5,fontsize=8)
    #ax.relim()
    #ax.autoscale()
    ax.set_xlim(-5,5)
    ax.figure.canvas.draw()
    ax = self.magplotWidget.ax
    ax.collections = []
    ax.lines = []
    L = []
    N = []
    if len(omag):
      L.append([ax.plot(mag[0],mag[1],label=mag[2],picker=5,marker='o',markerfacecolor='None',markersize=10,color='k',lw=0,zorder=11)[0] for mag in omag][0])
      N.append('Events')
    if len(tmag):
      L.append([ax.plot(mag[0],mag[1],label=mag[2],picker=5,marker='o',markerfacecolor='0.25',markersize=3,color='k',lw=0,zorder=11)[0] for mag in tmag][0])
      N.append('Stations (Used)')
    if len(tmagu):
      L.append([ax.plot(mag[0],mag[1],label=mag[2],picker=5,marker='o',markerfacecolor='None',markersize=1,color='k',lw=0,zorder=11)[0] for mag in tmagu][0])
      N.append('Stations (Unused)')
    ax.legend(L,N,frameon=False,loc=6,numpoints=1,fontsize=8)
    ax.plot([0,6],[0,6],'--',color='0.5',zorder=5)
    ax.relim()
    ax.autoscale()
    ax.figure.canvas.draw()

  def create_main_frame(self):
    'Create the main window widget with spliter, figure, toolbar'
    self.splitter = QSplitter(Qt.Horizontal) # main widget splitter
    self.treeSplitter = QSplitter(Qt.Vertical) # data view splitter
    self.splitter.addWidget(self.treeSplitter)
    # create a side by side layout for tables
    w = QWidget(parent=self.treeSplitter)
    self.tb = QHBoxLayout(w)
    self.tb.setAlignment(Qt.AlignLeft)
    self.treeSplitter.setStretchFactor(0,0)
    self.viewer = QWidget(self.splitter) # viewer for upper splitter panel
    self.vbox = QVBoxLayout(self.viewer) # vbox for figure and toolbar
    self.fig = Figure((6,3),dpi=100) # mpl figure (see matplotlib for details
    self.ax = self.fig.add_subplot(111) # mpl axes - see matplotlib for details
    self.ax.get_xaxis().get_major_formatter().set_useOffset(False) # make sure longitude are real numbers as %f
    self.ax.get_yaxis().get_major_formatter().set_useOffset(False) # make sure latitude are real numbers as %f
    self.ax.set_xlim(-160,90) # set starting view limits to map longitudes
    self.ax.set_ylim(-80,80) # set starting view limits to map latitudes
    self.ax.set_aspect('equal','datalim','C') # make sure lat and lon are not scaled differently and that widow size is fixed. data limits might change according to window size change
    self.ax.set_position([0,0,1,1]) # fill up the figure with the map axes
    self.ax.grid(True,color=[1,1,1,0.75]) # add grid
    self.ax.grid(GRIDON) # set grid to default state
    self.ax.tick_params('x',length=0,width=5,pad=-10,colors=[1,1,1,0.75]) # adjust x ticks
    self.ax.tick_params('y',length=0,width=5,pad=-20,colors=[1,1,1,0.75])# adjust y ticks
    [t.set_ha('left') for t in self.ax.xaxis.get_majorticklabels()] # adjust x ticks
    [t.set_va('bottom') for t in self.ax.yaxis.get_majorticklabels()] # adjust y ticks
    self.ax.drag_pan = self.drag_pan # replacing original drag_pan function of axes to home made one to make sure no panning out of geo bounderies
    self.canvas = FigureCanvas(self.fig) # set the canvas of figure
    self.canvas.setParent(self.viewer) # place canvas in viewer
    self.toolbar = NavigationToolbar(self.canvas, self.viewer) # add the toolbar of the canvas
    self.vbox.addWidget(self.canvas) # add canvas to layout
    self.vbox.addWidget(self.toolbar) # add toolbar to layout
    self.setCentralWidget(self.splitter) # set splitter as the main widget
    self.viewer.setFocus() # focus on viewer

  def create_menu(self):
    'Creates main menu'
    # Populate the menubar:
    ########### Add File submenu
    self.file_menu = self.menuBar().addMenu("&File")
    # load input data
    load_input_action = self.create_action("Load &E2 Log file",
            shortcut="Ctrl+O", slot=self.get_input_file,
            icon='document-open',tip="Load input data from a file(s)")
    # load reference data
    load_ref_action = self.create_action("Load &Reference",
            shortcut="Ctrl+R", slot=self.get_ref_file,
            icon='document-open',tip="Load reference data from a file(s)")
    # Quit
    quit_action = self.create_action("&Quit", slot=self.close,
            icon='system-shutdown',shortcut="Ctrl+Q", tip="Close the application")
    # populate the file submenu
    self.add_actions(self.file_menu,
            (load_input_action,load_ref_action, None, quit_action))
    ########### Add view submenue
    self.view_menu = self.menuBar().addMenu("&View")
    # toggle grid on or off
    togGrid_action = self.create_action("&Grid (on/off)",
            shortcut="Shift+Ctrl+G", slot=self.grid,
            icon=None,tip="Toggle map grid.",checkable=True)
    togGrid_action.setChecked(GRIDON)
    # populate the view submenu
    self.add_actions(self.view_menu,
            [togGrid_action])
    ########### Add help submenu
    self.help_menu = self.menuBar().addMenu("&Help")
    # Help
    help_action = self.create_action("&Help",
            shortcut='F1', slot=self.on_help,
            icon='Help',tip='help')
    # About
    about_action = self.create_action("&About",
            shortcut='F2', slot=self.on_about,
            tip='About This Application')
    # About QT
    aboutQt_action = self.create_action("&About QT",
            shortcut='F3', slot=self.on_aboutQt,
            tip='About QT')
    # License
    license_action = self.create_action("&License",
            shortcut='F4', slot=self.on_license,
            tip='Application License')
    # Populate help submenu
    self.add_actions(self.help_menu, (help_action,None,about_action,aboutQt_action,license_action))

  def updateMagLims(self):
    self.minmag = float(self.minmagLine.text())
    self.maxmag = float(self.maxmagLine.text())
    self.refreshbuttonHilight.setStrength(1)
    self.statusBar().showMessage('Magnitude limits updated',2000)

  def updatestarttime(self,starttime):
    self.starttime = starttime.toPyDateTime()
    self.refreshbuttonHilight.setStrength(1)
    self.statusBar().showMessage('Start time limit updated',2000)

  def updateendtime(self,endtime):
    self.endtime = endtime.toPyDateTime()
    self.refreshbuttonHilight.setStrength(1)
    self.statusBar().showMessage('End time limit updated',2000)

  def updateFixMags(self):
    self._fixmags = bool(self.fixmagsCheck.checkState())
    self.refreshbuttonHilight.setStrength(1)

  def updateFixDist(self):
    self._fixdist = bool(self.fixdistCheck.checkState())
    self.refreshbuttonHilight.setStrength(1)

  def updateReportedOnly(self):
    self.reportedOnly = bool(self.reportedOnlyCheck.checkState())
    self._associate = True
    self.refreshbuttonHilight.setStrength(1)

  def create_toolbox(self):
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    l = QGridLayout(w)
    self.refreshbutton = self.create_pushButton('Refresh', slot=self.refreshDisplay, shortcut=None, icon=None, tip='Update display after adjusting parameters')
    self.refreshbuttonHilight = QGraphicsColorizeEffect(self.refreshbutton)
    self.refreshbutton.setGraphicsEffect(self.refreshbuttonHilight)
    self.refreshbuttonHilight.setStrength(0)
    self.refreshbuttonHilight.setColor(QColor('orangered'))
    l.addWidget(self.refreshbutton,0,0)
    self.tb.addWidget(w)
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    l = QGridLayout(w)
    self.minmagLine = ConnectedLineEdit(str(self.minmag),0,self.maxmag)
    self.maxmagLine = ConnectedLineEdit(str(self.maxmag),self.minmag,10)
    self.minmagLine.setConnected(self.maxmagLine, 'bottom')
    self.maxmagLine.setConnected(self.minmagLine, 'top')
    l.addWidget(QLabel('Magnitude Limits'),0,1)
    l.addWidget(QLabel('Min'),1,0)
    l.addWidget(QLabel('Max'),2,0)
    l.addWidget(self.minmagLine,1,1)
    l.addWidget(self.maxmagLine,2,1)
    self.tb.addWidget(w)
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    v = QVBoxLayout(w)
    ww = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    v.addWidget(ww)
    l = QGridLayout(ww)
    W,E,S,N = self.bbox
    self.EAST = QLabel(str(E))
    self.WEST = QLabel(str(W))
    self.NORTH = QLabel(str(N))
    self.SOUTH = QLabel(str(S))
    self.bboxbutton = self.create_pushButton('Region', slot=self.Bbox.show, shortcut=None, icon=None, tip='Define a geographic region')
    l.addWidget(self.EAST,1,2,1,1,Qt.Alignment(Qt.AlignLeft))
    l.addWidget(self.WEST,1,0,1,1,Qt.Alignment(Qt.AlignRight))
    l.addWidget(self.NORTH,0,1,1,1,Qt.Alignment(Qt.AlignHCenter))
    l.addWidget(self.SOUTH,2,1,1,1,Qt.Alignment(Qt.AlignHCenter))
    v.addWidget(self.bboxbutton)
    self.tb.addWidget(w)
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    l = QGridLayout(w)
    self.starttimeLine = QDateTimeEdit()
    self.starttimeLine.setDateTime(QDateTime.fromMSecsSinceEpoch(self.starttime.timestamp*1000))
    self.starttimeLine.setDisplayFormat('yyyy-MM-dd HH:mm:ss')
    self.starttimeLine.setCalendarPopup(True)
    self.endtimeLine = QDateTimeEdit()
    self.endtimeLine.setDateTime(QDateTime.fromMSecsSinceEpoch(self.endtime.timestamp*1000))
    self.endtimeLine.setDisplayFormat('yyyy-MM-dd HH:mm:ss')
    self.endtimeLine.setCalendarPopup(True)
    l.addWidget(QLabel('Start:'),0,0)
    l.addWidget(self.starttimeLine,0,1)
    l.addWidget(QLabel('End:'),1,0)
    l.addWidget(self.endtimeLine,1,1)
    self.tb.addWidget(w)
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    l = QGridLayout(w)
    self.fixmagsCheck = QCheckBox('Fix Mags')
    self.fixmagsCheck.setChecked(self._fixmags)
    self.fixmagsCheck.setToolTip('Use Sadeh et al., (2014) Pd - magnitude relations to re-estimate magnitudes')
    l.addWidget(self.fixmagsCheck,0,0)
    self.tb.addWidget(w)
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    self.fixdistCheck = QCheckBox('Fix Dist')
    self.fixdistCheck.setChecked(self._fixdist)
    self.fixdistCheck.setToolTip('Use real distance to re-estimate magnitudes')
    l.addWidget(self.fixdistCheck,1,0)
    self.tb.addWidget(w)
    w = QWidget()
    w.setSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed)
    self.reportedOnlyCheck = QCheckBox('Reported Only')
    self.reportedOnlyCheck.setChecked(self.reportedOnly)
    self.reportedOnlyCheck.setToolTip('Use Only reported events from log file')
    l.addWidget(self.reportedOnlyCheck,2,0)

  def onBboxAccepted(self):
    '''get region limits from BboxForm.
       fires when Bbox is accepted'''
    if self.Bbox.validate(): # make sure limits are acceptable
      self.bbox = self.Bbox.getLims() # get limits
      w,e,s,n = self.bbox
      self.WEST.setText(str(w))
      self.EAST.setText(str(e))
      self.NORTH.setText(str(n))
      self.SOUTH.setText(str(s))
      self.refreshbuttonHilight.setStrength(1)

  def add_actions(self, target, actions):
    'Utility function for menu creation'
    for action in actions:
      if action is None:
        target.addSeparator()
      else:
        target.addAction(action)

  def create_action(  self, text, slot=None, shortcut=None,
                      icon=None, tip=None, checkable=False,
                      signal="triggered()"):
    'Utility function for menu actions creation'
    action = QAction(text, self)
    action.setIconVisibleInMenu(True)
    if icon is not None:
      i = QIcon.fromTheme(icon,QIcon(":/%s.png" % icon))
      action.setIcon(i)
    if shortcut is not None:
      action.setShortcut(shortcut)
    if tip is not None:
      action.setToolTip(tip)
      action.setStatusTip(tip)
    if slot is not None:
      self.connect(action, SIGNAL(signal), slot)
    if checkable:
      action.setCheckable(True)
    return action

  def create_pushButton(self,text,toolbar=None, slot=None, shortcut=None, icon=None, tip=None):
    'Utility function for button creation'
    # create the button
    button = QPushButton(text,self)
    # populate properties
    if slot:
      # connect a function
      button.clicked.connect(slot)
    if icon:
      # add icon
      i = QIcon.fromTheme(icon,QIcon(":/%s.png" % icon))
      button.setIcon(i)
      button.setIconSize(QSize(24,24))
    if shortcut:
      # set the shortcut
      button.setShortcut(shortcut)
    if tip:
      # add tooltip and status tip
      button.setToolTip(tip)
      button.setStatusTip(tip)
    if toolbar:
      # add the button to a toolbar (or any widget)
      toolbar.addWidget(button)
    return button

  def create_status_bar(self):
    'Add a status bar'
    # set default message
    self.status_text = QLabel("Ready")
    self.connstatLabel = QLabel()
    self.statusBar().addWidget(self.status_text, 1)
    self.statusBar().addPermanentWidget(self.connstatLabel)


  def on_about(self):
    'show a messagebox about the application'
    msg = "<p align='center'><big>E2ReviewTool</big><br><br> \
    Review ElarmS E2 Results<br><br> \
    <small>Created<br> \
    by<br> \
    Ran Novitsky Nof @ BSL, 2015</small><br><br>\
    <a href='http://ran.rnof.info/'>http://ran.rnof.info</a><p>"
    QMessageBox.about(self,"About", msg.strip())

  def on_aboutQt(self):
    'show a messagebox about QT'
    QMessageBox.aboutQt(self,'')

  def on_license(self):
    'GPL licanse message'
    msg = "<p><b>This</b> is a free software; you can redistribute it and/or modify it under the \
terms of the GNU General Public License as published by the Free Software \
Foundation; either version 3 of the License, or (at your option) any later \
version.</p>\
<p><b>This application</b> is distributed in the hope that it will be useful, but WITHOUT ANY \
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR \
A PARTICULAR PURPOSE.  See the GNU General Public License for more details.</p> \
<p>You should have received a copy of the GNU General Public License along with \
this application; if not, see <a href='http://www.gnu.org/licenses/'>http://www.gnu.org/licenses/</a>.</p>"
    QMessageBox.about(self,"Application Licanse", msg.strip())

  def on_help(self):
    'Show help on a message window. Uses the argparse help'
    msg = '<pre>'+parser.format_help()#.replace('\n','</p><p>')
    msg = msg.replace('ran.nof@gmail.com',"<a href='mailto:ran.nof@gmail.com'>ran.nof@gmail.com</a>").strip()+'<\pre>'
    QMessageBox.about(self,"Help", msg)

  def message(self,msg,title='Error'):
    'a simple message window'
    QMessageBox.about(self,title,msg)

def main(args):
  # create the application
  app = QApplication(sys.argv)
  #splash
  splash_pix = QPixmap('splash.png')
  splash_pix.scaledToHeight(50)
  splash = QSplashScreen(splash_pix)
  splash.setMask(splash_pix.mask())
  splash.show()
  app.processEvents()
  # populate the QT4 form
  appwin = AppForm(splash,args)
  appwin.show()
  splash.finish(appwin)
  # run the application
  sys.exit(app.exec_())

if __name__=="__main__":
  # parse the arguments
  args = parser.parse_args(sys.argv[1:])
  mpl.rcParams['font.size']=float(FONTSIZE)
  mpl.rcParams['figure.facecolor']='w'
  main(args)
else:
  print '''

import glob
app = QApplication([])
ifiles = ifiles = glob.glob('/work/work/GIIE2LOGS/E2_20150[567]*.log')+['/work/work/E2r_20151029_201201-201504.log']
rfiles = ['/work/work/GII_EQ_20150416_20150801.csv','/work/work/GII_EQ_201201-201504.csv']
args = parser.parse_args(['-i']+ifiles+['-r']+rfiles)
self = AppForm(None,args)
self.show()



  '''


w='''
d1 = []
d2 = []
d3 = []
d4 = []
for EID in edict:
  event = edict[EID]
  refID = event['refID']
  origins = [origin for origin in odict[refID] if origin['first']]
  if not len(origins): continue
  origin = origins[0]
  triggers = tdict['.'.join([refID,str(origin['ver'])])]
  fig = figure(1)
  ax1 = fig.add_subplot(221)
  ax1.set_title('ElarmS')
  ax2 = fig.add_subplot(222)
  ax2.set_title('Corrected distance')
  ax3 = fig.add_subplot(223)
  ax3.set_title('Sadeh et al., 2013 (dist. corr.)')
  ax4 = fig.add_subplot(224)
  ax4.set_title('Sadeh et al., 2013')
  for trigger in triggers:
    m0 = event['mag']
    m1 = trigger['pdmag']
    m2 = pdmag(trigger['logpd'],geo.inv(event['lon'],event['lat'],trigger['lon'],trigger['lat'])[-1]/1000.)
    m3 = pdmag1(trigger['logpd'],geo.inv(event['lon'],event['lat'],trigger['lon'],trigger['lat'])[-1]/1000.)
    m4 = pdmag1(trigger['logpd'],trigger['distkm'])
    print trigger['sta']+'.'+trigger['chn'],m0,m1,m2,m3,m4,m0-m1,m0-m2,m0-m3,m0-m4
    if not 'MMA0' in trigger['sta']:
      ax1.scatter(m0,m1)
      ax2.scatter(m0,m2)
      ax3.scatter(m0,m3)
      ax4.scatter(m0,m4)
      d1.append(m0-m1)
      d2.append(m0-m2)
      d3.append(m0-m3)
      d4.append(m0-m4)
'''


def pdtp(logpds,logtps):
  return log10(mean([10**i for i in logpds])),log10(mean([10**i for i in logtps]))

def on_point_Pick(event):
  try:
    event.canvas.toolbar.set_message(event.artist.get_label())
  except Exception as E:
    print str(E)

'''
fig = figure()
ax = gca()
fig.canvas.mpl_connect('pick_event',on_point_Pick)
for i in xrange(self.eventWidget.topLevelItemCount()):
  item = self.eventWidget.topLevelItem(i)
  if not item.preffered: continue
  origin = item.preffered
#  for tr in origin.triggers:
#    if tr.updm and 'PRNI' in tr.ID:
#      tr.logpd=log10(10**tr.logpd)/2.)
  if not any([tr.dist<=100 for tr in origin.triggers]): continue
  logpds = array([tr.logpd for tr in origin.triggers if tr.updm and not 'MBRI' in tr.ID and not 'DAM2' in tr.ID])
  pdsnrs = array([tr.pdsnr for tr in origin.triggers if tr.updm and not 'MBRI' in tr.ID and not 'DAM2' in tr.ID])
  logtps = array([tr.logtp for tr in origin.triggers if tr.utpm])
  tpsnrs = array([tr.tpsnr for tr in origin.triggers if tr.utpm])
  logpd = log10((10**logpds*pdsnrs).sum()/pdsnrs.sum())
  logtp = log10((10**logtps*tpsnrs).sum()/tpsnrs.sum())
  if sum(pdsnrs>100)<1:
    plot(logtp,logpd,'.',c=['b','r','r','b','k'][origin.event.typ],label=origin.EID,picker=5,alpha=0.2)
  else:
    plot(logtp,logpd,'.',c=['b','r','r','b','k'][origin.event.typ],label=origin.EID,picker=5)
x=arange(-1,1,0.2)
plot(x,-3.728+x*2.817,'k--')
fig.canvas.mpl_connect('pick_event',on_point_Pick)
'''




























