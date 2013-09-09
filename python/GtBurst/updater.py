#!/usr/bin/env python

import urllib, os
from GtBurst.GtBurstException import GtBurstException
import GtBurst
import md5

#Set a global timeout of 10 seconds for all web connections
import socket
socket.setdefaulttimeout(60)

def update(debug=False):
  #Download file_list file
  try:
    os.remove("__file_list")
  except:
    pass
  urllib.urlcleanup()
  try:
    urllib.urlretrieve("http://sourceforge.net/projects/gtburst/files/Stable/__file_list/download","__file_list")
  except socket.timeout:
    raise GtBurstException(11,"Time out when connecting to Sourceforge. Check your internet connection, or that you can access http://sourceforge.net/projects/gtburst/files/, then retry")
  except:
    raise GtBurstException(1,"Problems with the download. Check your connection, and that you can reach http://sourceforge.net/projects/gtburst/files/, then retry")
  pass
  
  f                           = open('__file_list')
  files                       = f.readlines()
  f.close()
  os.remove("__file_list")
  
  #Get the path of the gtburst installation
  path                      = GtBurst.__file__
  installationPath          = os.path.join(os.path.sep.join(path.split(os.path.sep)[0:-3]))
  
  nUpdates                  = 0
  for ff in files:
    atoms                     = ff.split()
    pathname                  = atoms[-1].replace('*','')
    if(ff[0]=='d'):
      #This is a directory, skipping
      if(debug):
        print("Skip directory %s" %(pathname))
    elif(ff.find("__file_list")>=0):
      if(debug):
        print("Skipping %s..." %(ff))
    else:
      remoteMD5                    = atoms[0]
      if(debug):
        print("File %s has remote MD5 checksum %s" %(pathname,remoteMD5))
      
      #Get the MD5 of the same file in the GtBurst package path
      pathnameThisSys         = pathname.replace("/",os.path.sep)
      localPath               = os.path.join(installationPath,pathnameThisSys)
      if(not os.path.exists(localPath)):
        print("File %s does not exist in the current installation. Creating it..." %(localPath))
        downloadFile(pathname,localPath)
        nUpdates              += 1
      else:
        #File exists. Check the size
        localMD5             = md5.md5(open(localPath, 'rb').read()).hexdigest()
        if(localMD5!=remoteMD5):
          print("Updating %s..." %(localPath))
          downloadFile(pathname,localPath)
          nUpdates              += 1
        else:
          if(debug):
            print("NOT updating %s (local MD5: %s, remote MD5: %s)..." %(localPath,localMD5,remoteMD5))
          pass
        pass
        if(debug):
          print("\n\n")
    pass
  pass
  
  return nUpdates
pass

def downloadFile(remotepathname,localpathname):
  try:
    urllib.urlretrieve("http://sourceforge.net/projects/gtburst/files/Stable/%s/download" % remotepathname,localpathname)
  except socket.timeout:
    raise GtBurstException(11,"Time out when connecting to Sourceforge. Check your internet connection, or that you can access http://sourceforge.net/projects/gtburst/files/, then retry")
  except:
    raise GtBurstException(1,"Problems with the download. Check your connection, and that you can reach http://sourceforge.net/projects/gtburst/files/, then retry")
  pass
  
