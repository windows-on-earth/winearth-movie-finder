## -----------------------------------------------------------------------------
## Copyright(c) 2010 - 2024 ViSUS L.L.C.,
## Scientific Computing and Imaging Institute of the University of Utah
##
## ViSUS L.L.C., 50 W.Broadway, Ste. 300, 84101 - 2044 Salt Lake City, UT
## University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
##
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met :
##
## * Redistributions of source code must retain the above copyright notice, this
## list of conditions and the following disclaimer.
##
## * Redistributions in binary form must reproduce the above copyright notice,
## this list of conditions and the following disclaimer in the documentation
## and/or other materials provided with the distribution.
##
## * Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## For additional information about this project contact : pascucci@acm.org
## For support : amy@visus.net
## -----------------------------------------------------------------------------

## This software opens a GUI for evaluating the contrast stretching of photographs
## The algorithm is as follows:
## Upload the image,
## Downsample the image
## compute the histogram,
## remove an offset from the histogram (emphrically chosen as 75, may need to alter this if we do histogram on full res)
## find the edges of the new histogram
## Use those values for the new contrast range on the full size image.
##
## There is an option to display small images on a grid for seeing what's happening
##
## And an option for just applying the algorithm
##
## Then the GUI class displays two images side by side, first the original, next the altered images
## The two sliders below indicate standard deviations on each side of the mean
##

## defaults.USE_PREVIEW_MODE  = True
## resizes the images to 300x200 (with aspect ratio preserved) in order to save on computation.
## Commandline version (which is coming soon) will use the full resolution


import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import traceback

import scipy.stats as sts

from PIL import Image
from PIL import ImageCms

import panel as pn
pn.extension()
import numpy as np
from PIL import Image
from bokeh.io import curdoc
from matplotlib.figure import Figure
import json
import argparse


class CLIApp:
    def __init__(self):
        super().__init__()
        self.fileOptions = [];
        self.curFile = None
        self.curDir = None
        self.outDir = None
        self.sourceImgFileName = None
        self.dstImgFileName = None

        self.srcImageObject = None
        self.srcImageObjectSmall = None
        self.dstImageObject = None
        self.dstImageObjectSmall = None

        self.minI= 86
        self.maxI= 230
        self.minO    = 0
        self.maxO    = 255
        self.minBin0 = 0
        self.maxBin0 = 255
        self.max_bin = 128
        self.amean, self.avar, self.astdev = 128,245,15.0
        self.myData =None

        self.fig = Figure(figsize=(3, 1.5))

    def setData(self):
        self.myData = {
            'filename':self.curFile,
            'minI': self.minI,
            'maxI': self.maxI,
            'minbin': self.minBin0,
            'maxbin': self.maxBin0,
            'path': self.curDir,
            'astdev':self.astdev,
            'amean':self.amean,
        }


    def find_transition_edges(self, n, bins):
        # print(*n, sep = ',')
        n -= 75
        n[n < 0] = 0
        left_edge = 0
        right_edge = 255
        # print(*n, sep = ',')
        # Find the first edge where the count goes from 0 to non-zero
        for i in range(0, len(n) - 1, 1):
            if n[i]>100:
                left_edge = 0
                break
            elif n[i] < 0.01 and n[i + 1] >= 1:
                # print('found left at {} {}'.format( n[i], n[i+1]))
                left_edge = bins[i + 1]
                break

        # Find the first edge where the count goes from non-zero to 0
        for i in range(len(n)-2, 0, -1):
            if n[len(n)-1] > 100:
                right_edge = len(n)-1
                break
            elif n[i+1 ] < 0.01 and n[i] >= 1:
                # print('found right at {} {}'.format( n[i+1], n[i]))
                right_edge = bins[i+1]
                break
        # print('Edges: {} to {}'.format(left_edge, right_edge))
        return int(left_edge), int(right_edge)


    def calcMeasures(self):
        a = np.array(self.srcImageObjectSmall.getdata())

        n, bins, patches = plt.hist(a.ravel(),255,[0,254])
        mids = 0.5*(bins[1:] + bins[:-1])
        self.amean = np.average(mids, weights=n)
        self.avar = np.average((mids - self.amean)**2, weights=n)
        self.astdev =np.sqrt(self.avar)

        self.minBin0 = 0
        self.maxBin0 = 255
        # Find the edges
        left_edge, right_edge = self.find_transition_edges(n, bins)
        #print("Find Transition:  left_edge={} right edge{}".format(left_edge,right_edge))

        if (left_edge and left_edge >= 0 and left_edge <= 255):
            self.minBin0 = left_edge
        if (right_edge and right_edge >= 0 and right_edge <= 255):
            self.maxBin0 = right_edge
        if (self.minBin0 > self.maxBin0):
            self.minBin0 = 0
        if (self.maxBin0 < self.minBin0):
            self.maxBin0 = 255

        self.setData()
        # print("mean")
        # print(self.amean)
        # print("var")
        # print(self.avar)
        # print("stdev")
        # print(self.astdev)

        #return (self.amean, self.avar, self.astdev)
    def normalizeChannelPointOp(self, intensity):
        iI      = intensity
        d = self.maxI-self.minI
        if (d== 0):
            d = .001
        iO = (iI-self.minI)*(((self.maxO-self.minO)/(d))+self.minO)
        return iO


    def applyContrastStretchingSrcDst(self, src  ):
        #print('apply Contrast Stretching')
        # Split the red, green and blue bands from the Image

        multiBands      =  src.split()

        # Apply point operations that does contrast stretching on each color band
        normalizedRedBand      = multiBands[0].point(self.normalizeChannelPointOp)
        normalizedGreenBand    = multiBands[1].point(self.normalizeChannelPointOp)
        normalizedBlueBand     = multiBands[2].point(self.normalizeChannelPointOp)

        # Create a new image from the contrast stretched red, green and blue brands
        dstImageObject  = Image.merge("RGB", (normalizedRedBand, normalizedGreenBand, normalizedBlueBand))


        return dstImageObject
    def applyContrastStretching(self ):
        self.dstImageObject = self.applyContrastStretchingSrcDst( self.srcImageObject )
        self.dstImageObjectSmall = self.applyContrastStretchingSrcDst(self.srcImageObjectSmall )

    def getPlot(self):
        self.fig.clf()
        #print('plot')

        # adding the subplot
        plot1 = self.fig.add_subplot(111)
        # plotting the graph
        plot1.vlines(self.minI, -10, 5000, colors='black', linestyles='-', label='min', data=None )
        plot1.vlines(self.maxI, -10, 5000, colors='black', linestyles='-', label='min',  data=None )

        a = np.array(self.srcImageObjectSmall.getdata())
        n, b, patches = plot1.hist(a.ravel(),255,[0,254])
        plot1.plot()
        plot1.axes.xaxis.set_ticklabels([])
        plot1.axes.yaxis.set_ticklabels([])

        plot1.text(10, 7020, 'minI: {}= n {}'.format(self.minI, n[self.minI]))
        plot1.text(10, 6020, 'maxI: {}= n {}'.format(self.maxI, n[self.maxI]))


    def resize_imageObjects(self, imageObject):
        if (imageObject):
            imageObject = imageObject.resize((300, 200), Image.LANCZOS);
        return imageObject

    def open_image(self, filename, filepath, outdir):
        try:
            self.curFile = filename
            self.curDir = filepath
            self.outDir = outdir
            self.minI= 86
            self.maxI= 230
            self.srcImageObject     =None

            self.dstImageObject     = None
            self.sourceImgFileName =   os.path.join(self.curDir ,self.curFile)
            self.dstImgFileName =   os.path.join(self.curDir ,'contrast_'+self.curFile)
            self.srcImageObject = Image.open(self.sourceImgFileName)
            self.srcImageObjectSmall = self.resize_imageObjects(self.srcImageObject )
            self.calcMeasures()
            self.minI= self.minBin0 #+ int(round(self.max_contrast_val))
            self.maxI= self.maxBin0 #+ int(round( self.min_contrast_val))
            self.setData()
            self.applyContrastStretching()
        except IOError as e:
            print ("I/O error({0}): {1}".format(e.errno, e.strerror))
        except: #handle other exceptions such as attribute errors
            print ("Unexpected error:", sys.exc_info()[0])
    def save(self):
        newFilename = '{}/contrast_{}.png'.format( self.outDir, self.curFile)
        self.dstImageObject.save(newFilename)
        return newFilename


#### Get arguments
# Initialize parser
parser = argparse.ArgumentParser(description='Contrast stretching for Windows on the Earth ISS Photograph ',
                                 epilog='usage: python contrastForWOE.py -p <src img path> \n optional: \n  -o <output dir> \n  -f filename to process just one file')
parser.add_argument('-p', '--path',    type=str,
                    help='path to image file, not including image name', required=True, default=
                    'images/')
parser.add_argument('-f', '--filename',    type=str,
                    help='if a file is given, then it will only processs the file; image file name (with extension, no path)', required=False)
parser.add_argument('-o', '--output',    type=str,
                    help='output directory, if not provided it will be stored in same directory as source images', required=False, default= 'save')
parser.add_argument('-v', '--viewweb',   action="store_true",
                    help='view web page of images', required=False )
parser.add_argument('-d', '--donotsave',   action="store_true",
                    help='do not save images, just show web page', required=False )


def processFile(filename, originPath,output ):
    if filename.lower().endswith(('.jpg', '.png', '.tif', '.tiff')):
        app = CLIApp()
        app.curFile = filename
        app.curDir = originPath
        app.outDir =  output
        app.open_image(app.curFile,  app.curDir, app.outDir )
        outputSaved = app.save()
        return outputSaved



args = parser.parse_args()
originPath = args.path
filename=None
outputpath=None
fullpathfile = None

try:
    if args.filename is not None:
        filename = args.filename
        fullpathfile  = os.path.join(originPath, filename)

    if args.output is not None:
        outputpath = args.output
    else:
        outputpath = args.path
except:
    print ("Unexpected error:", sys.exc_info()[0])
    print('ERROR: processing inputs')

if ( fullpathfile and os.path.isfile(fullpathfile) and os.access(fullpathfile, os.R_OK)):
    #Parse only a file
    try:
        o = processFile(filename, originPath, args.output)
        print('processed file and saved to {}'.format(o))
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print('ERROR: Failed to process image: {}/{}'.format(originPath,x))

else:
    #Parse a whole directory of files
    print('Processing whole directory: {}'.format(originPath))
    fileOptions = []
    gspec = pn.GridSpec()
    i = 0
    for x in os.listdir(originPath):
        print("processing {}/{}".format(originPath,x))
        try:
            if x.lower().endswith(('.jpg', '.png', '.tif', '.tiff')):
                app = CLIApp()
                app.curFile = x
                app.curDir = originPath
                app.outDir = args.output
                app.open_image(app.curFile,  app.curDir, app.outDir )
                if (args.viewweb and gspec):
                    app.getPlot()
                    gspec[i,   0  ] =pn.pane.Image(app.srcImageObjectSmall)
                    gspec[i,   1  ] =pn.pane.Image(app.dstImageObjectSmall)
                    gspec[i,   2  ] =pn.pane.Matplotlib(app.fig, dpi=144)
                    gspec[i,   3  ] =pn.pane.JSON(app.myData, name='Current Settings')
                if not args.donotsave:
                    o = app.save()
                    print('processed file and saved to {}'.format(o))
                fileOptions.append(x)
        except IOError as e:
            print ("I/O error({0}): {1}".format(e.errno, e.strerror))
            print('ERROR: Failed to process image: {}/{}'.format(originPath,x))
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
            print(traceback.format_exc())
            print('ERROR: Failed to process image: {}/{}'.format(originPath,x))
        except: #handle other exceptions such as attribute errors
            print ("Unexpected error:", sys.exc_info()[0])
            print('ERROR: Failed to process image: {}/{}'.format(originPath,x))
        i+=1
    if (args.viewweb):
        pn.serve(gspec)