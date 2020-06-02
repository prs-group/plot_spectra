#!/usr/bin/env python3
# vim:tabstop=4:softtabstop=4:shiftwidth=4

import sys
from argparse import ArgumentParser
from argparse import FileType
from argparse import Action
import argparse
from plot_functions import *
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib import widgets
import matplotlib.image as image
import re

rc('text', usetex=True)

# Some tweaks
topValueForNumberOfLabels={ "0": 1.0, 
        "1": 0.94, 
        "2": 0.89,
        "3" : 0.85,
        "4" : 0.8 }

class AppendRange(Action):
    def __init__(self,option_strings,dest,nargs='*',**kwargs):
        super(AppendRange,self).__init__(option_strings,dest,nargs,**kwargs)
    
    def __call__(self,parser,namespace,values,option_string=None):
        rangeList=[]
        rangeRegex=re.compile(r'([+\-]?\d+\.?\d+)-([+\-]?\d+\.?\d+)')
        itValues=values
        if any(isinstance(el,list) for el in values): 
            itValues=values[0]

        for value in itValues:
            result=rangeRegex.search(value)
            if result:
                rangeList.append((float(result.group(1)),float(result.group(2))))
            else:
                raise(ValueError("You did not specify a valid range"))


        items = getattr(namespace,self.dest,None)
        if items is None:
            items = []
            
        for item in rangeList:
            items.append(item)

        setattr(namespace,self.dest,items)

def RangesY(value):
    values=value.split()
    if len(values) != 2:
        raise argparse.ArgumentError()

    return values

def Ranges(value):
    values=value.split()
    return values

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def ListOfFloats(v):
    if any(isinstance(el,float) for el in v):
        return v
    
    floats=v.split()
    listOfFloats=[]
    for f in floats:

        if f.lower() in ('-','nan','none','notfound'):
            listOfFloats.append(None)

        else:
            listOfFloats.append(float(f))
    
    return listOfFloats



def _generateParser():
    parser=ArgumentParser(description='Plot an IR spectra.')
    plotOptionGroup=parser.add_argument_group(title='Plot Options',
            description='Options that affect the whole plot')
    spectraOptionGroup=parser.add_argument_group(
            title='Options that affect the each spectrum',
            description="""These options will influence the settings for each graph. 
            If not mentioned explicitly each option can be specified multiple times
            for each spectrum""")
    fileOptionGroup=parser.add_argument_group(title='Spectra and Output File',
            description="Specifiy as many spectra. The last argument is the output figure")

    fileOptionGroup.add_argument('spectraDataFiles',nargs='+',type=FileType('r'),
            action='store')
    fileOptionGroup.add_argument('outputFile',nargs=1,type=FileType('w+'))

    plotOptionGroup.add_argument('--distanceSubplots',type=float,nargs='*',
            default=0.03,action="append",
            help="The padding between the axis breaks.")
    plotOptionGroup.add_argument('--labelXAxis',type=str,nargs=1,
            action="store",
            default=r"wavenumber [cm$^{-1}$]",
            help="The label for the X Axis")
    plotOptionGroup.add_argument('--labelYAxis',type=str,nargs=2,
            default=[r'absorption [arb. units]',r'intensity [km\ mol$^{-1}$]'],
            action="append",
            help="The label for the Y Axis")
    plotOptionGroup.add_argument('--majorTicksX',type=float,nargs='*',
            action="store",
            help="Set the minor ticks on the x axis.")
    plotOptionGroup.add_argument('--minorTicksX',type=float,nargs='*',
            action="store",
            help="Set the minor ticks on the x axis.")
    plotOptionGroup.add_argument('--plotLimitsX',nargs=1,
            action=AppendRange,type=Ranges,
            help="Set the range of the x axis by specifying the range (e.g. 4000-2000)")
    plotOptionGroup.add_argument('--interactive',
            action="store_true",
            help="Set the range of the x axis by specifying the range (e.g. 4000-2000)")
    plotOptionGroup.add_argument('--plotLimitsY',nargs=1,
            action=AppendRange,type=RangesY,
            help="Set the range of the y axis by specifying the range (e.g. 4000-2000)")
    plotOptionGroup.add_argument('--image', nargs=1,
            type=FileType('r'),action='store',help="Show image in plot")
    plotOptionGroup.add_argument('--imagePosition',type=float,nargs=4,
            action="store", default=[0.0,0.0,0.2,0.2],
            help="4 floats: left,bottom,width,height. Starting from north west.")

    spectraOptionGroup.add_argument('--colors',nargs='*',
            action='store',
            help="Set the color for each plot",default="k")
    spectraOptionGroup.add_argument('--labels',nargs='*',
            action='store',
            help='The labels for each plot')
    spectraOptionGroup.add_argument('--assignments',nargs='*',
            type=ListOfFloats,action='store',
            help='If specified for two graphs assignment lines are drawn')
    spectraOptionGroup.add_argument('--invert',action="store",
            default="false",nargs='*',type=str2bool,
            help="Specify if plot should be inverted")
    spectraOptionGroup.add_argument('--yshift',type=float,
            nargs='*',action="store",default=0.0,
            help="Shift the data along the y axis")
    spectraOptionGroup.add_argument('--peakWidth',type=float,
            nargs='*',action="store",default=5.0,
            help="The width of the peaks for assignment.")
    spectraOptionGroup.add_argument('--yaxisIndex', type=int,
            nargs='*', action="store",default=0,choices=[0,1],
            help="Select the y axis on which to plot.")

    return parser

def _getValueForPlotOptionAtPosition(option,position):

    if not isinstance(option,list):
        return option

    if len(option) < position:
        print("Warning: option for index {} could not be found!".format(position))
        return option[0]
    else:
        return option[position]

def _compileSpectraOptionsFromOptions(parsedOptions,file,index):
    plotOptions={}
    # If the option for index does not exist take the first one
    # as it is the default option
    colors=(_getValueForPlotOptionAtPosition(parsedOptions.colors,index),)
    labels=(_getValueForPlotOptionAtPosition(parsedOptions.labels,index),)
    assignments=(_getValueForPlotOptionAtPosition(parsedOptions.assignments,index),)
    invert=(_getValueForPlotOptionAtPosition(parsedOptions.invert,index),)
    yshift=(_getValueForPlotOptionAtPosition(parsedOptions.yshift,index),)
    peakWidth=_getValueForPlotOptionAtPosition(parsedOptions.peakWidth,index)

    plotOptions["filenames"]=(file,)
    plotOptions["colors"]=colors
    plotOptions["labels"]=labels
    plotOptions["assignments"]=assignments
    plotOptions["invert"]=invert[0]
    plotOptions["yshift"]=yshift
    plotOptions["peakWidth"]=peakWidth

    return plotOptions

def _compilePlotOptions(parsedOptions):
    plotOptions={}
    plotOptions["distanceBetweenPlots"]=parsedOptions.distanceSubplots
    plotOptions["labelsXAxis"]=(parsedOptions.labelXAxis,)
    plotOptions["labelsYAxis"]=parsedOptions.labelYAxis
    plotOptions["majorTicksX"]=parsedOptions.majorTicksX
    plotOptions["minorTicksX"]=parsedOptions.minorTicksX
    plotOptions["plotLimitsX"]=parsedOptions.plotLimitsX
    plotOptions["plotLimitsY"]=parsedOptions.plotLimitsY

    return plotOptions

def main():
    parser=_generateParser()
    arguments=parser.parse_args()
    plotOptions=_compilePlotOptions(arguments)
    # The number of axis interruptions is determined by length of plotLimitsX
    numberOfSubplots=1
    if plotOptions["plotLimitsX"]:
        numberOfSubplots=len(plotOptions["plotLimitsX"])

    figure,axes=setupPlot(numberOfSubplots,plt,**plotOptions,tight=False)
    plt.subplots_adjust(bottom=0.1,left=0.2)

    plots=[]
    assignmentPoints=[]

    for index,file in enumerate(arguments.spectraDataFiles):
        spectraOption=_compileSpectraOptionsFromOptions(arguments,file,index)
        if type(arguments.yaxisIndex) == list:
            yaxes=axes[arguments.yaxisIndex[index]]
        else:
            yaxes=axes[arguments.yaxisIndex]

        theplot,assignmentPoint=plot(**spectraOption,figure=figure,
                axes=yaxes)
        plots.append(theplot[0])
        assignmentPoints.append(assignmentPoint)

    # By default we merge the last groups as they are computational bands
    if assignmentPoints:
        experimentalBands=assignmentPoints[0]
        computedBands=[]
        for bands in assignmentPoints[1:]:
            computedBands+=bands
        
        if len(experimentalBands) != len(computedBands):
            print("Warning: Experimental and computation band length do not match")

        drawAssignments(experimentalBands,computedBands,axes[1])



    if arguments.interactive:
        plt.subplots_adjust(bottom=0.25,left=0.2)
        axes[0][0].margins(x=1.0)
        axisY1Slider=figure.add_axes([0.15,0.15,0.65,0.03])
        axisY2MaxSlider=figure.add_axes([0.15,0.19,0.30,0.03])
        axisY2MinSlider=figure.add_axes([0.6,0.19,0.30,0.03])
        y2Min,y2Max=axes[1][0].get_ylim()

        margin=0.2
        y2MinUpper=y2Min+abs(y2Min*margin)
        y2MinLower=y2Min-abs(y2Min*margin)
        y2MaxUpper=y2Max+abs(y2Max*margin)
        y2MaxLower=y2Max-abs(y2Max*margin)

        step=int(y2Min/100)

        y1slider=widgets.Slider(ax=axisY1Slider,label="Y1 Margin",
            valmin=-0.5,valmax=5.0,valinit=1.0,
            valstep=0.01,orientation='horizontal')

        y2MaxSlider=widgets.Slider(ax=axisY2MaxSlider,label="Y2 Max",
            valmin=y2MaxLower,valmax=y2MaxUpper,valinit=y2Max,
            valstep=step,orientation='horizontal')

        y2MinSlider=widgets.Slider(ax=axisY2MinSlider,label="Y2 Min",
            valmin=y2MinLower,valmax=y2MinUpper,valinit=y2Min,
            valstep=step,orientation='horizontal')

        def updateY1(val):
            maxY = y1slider.val
            for axis in axes[0]:
                axis.margins(y=maxY)
            figure.canvas.draw_idle()

        def updateY2(val):
            maxY = y2MaxSlider.val
            minY = y2MinSlider.val
            for axis in axes[1]:
                axis.set_ylim((minY,maxY))

            figure.canvas.draw_idle()

        y1slider.on_changed(updateY1)
        y2MaxSlider.on_changed(updateY2)
        y2MinSlider.on_changed(updateY2)
    

    labels= [ l.get_label() for l in plots ]

    # Insert Picture
    if arguments.image:
        insetPictureArray=image.imread(arguments.image[0].name)
        newax=figure.add_axes(arguments.imagePosition,anchor='NW',zorder=5)
        newax.imshow(insetPictureArray)
        newax.axis('Off')


    plt.figlegend(plots,labels,frameon=False,loc="upper right")
    plt.subplots_adjust(bottom=0.1,
            top=topValueForNumberOfLabels[str(len(labels))],
            left=0.12,
            right=0.88,
            wspace=0.05)
    #plt.show()
    figure.savefig(arguments.outputFile[0].name,dpi=1200)
        

if __name__ == "__main__":
    main()

