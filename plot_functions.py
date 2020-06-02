#!/usr/bin/env python3
# vim:tabstop=4:softtabstop=4:shiftwidth=4

import matplotlib.lines as lines
import matplotlib.ticker as ticker
from scipy.signal import find_peaks
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import math
import pathlib

_suffixToDelimiter={ ".csv" : ";", ".tsv": " ", ".txt": " ", ".dpt" : "\t" }

def getPeakYForX(xdata,ydata,x,offset=1.5):
    # First sort the x and ydata
    xydata=[ [x,y] for x,y in zip(xdata,ydata) ]
    xydata.sort(key=lambda point: point[0])
    # Take a subset of the data
    lowerLimit=x-offset
    upperLimit=x+offset
    condlist= (xdata>=lowerLimit) & (xdata<=upperLimit)
    selectedYValues=np.extract(condlist,ydata)
    indexHigh=find_peaks(selectedYValues)[0]
    indexLow=find_peaks(-1*selectedYValues)[0]

    if indexLow.size > 0:
        return selectedYValues[indexLow[0]],True
    elif indexHigh.size > 0:
        return selectedYValues[indexHigh[0]],False
    else:
        return None,False

def importData(filename,delim=None):

    if delim:
        xdata,ydata=np.loadtxt(filename,delimiter=delim,unpack=True)
    else:
        # For a filetype 
        suffix=""
        if hasattr(filename,"name"):
            suffix=pathlib.Path(filename.name).suffix
        elif hasattr(filename,"suffix"):
            suffix=filename.suffix
            # If it is a string
        else:
            suffix=pathlib.Path(filename).suffix

        # Map the delimiter of the file
        try:
            xdata,ydata=np.loadtxt(filename,delimiter=_suffixToDelimiter[suffix],unpack=True)
        except ValueError as e:
            print("The file {} could not be read!".format(filename))
            raise e
        except KeyError:
            raise KeyError("The suffix {} is not known! Specifiy suffix explicitly".format(suffix))
        except ValueError as e:
            raise e


    return xdata,ydata

def removeSpines(ax):
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

def _createFigureAndAxes(numberOfXAxisInterruptions,tight=True):
    """
    Creates the subplots and axes for an x-axis interrupted plot.
    For each subplot two y-axis are provided and return seperatly.

    Returns: The figure,list of Y1Axes, list of Y2Axes
    """
    # For each interruption the graph is expanded to the right
    axes1=[]
    figure,axes1=plt.subplots(1,numberOfXAxisInterruptions,sharey=True)
    if numberOfXAxisInterruptions == 1:
        axes1=[ axes1 ]
    figure.set_tight_layout(tight)
    # We setup the second y-axis
    axes2=[ axis.twinx() for axis in axes1 ]
    return figure,axes1,axes2

def _setAxesRanges(axes1,axes2,plotLimitsX,plotLimitsY):
    """
    Sets the Axis Limits for X and Y axis using a list of axes
    """
    # The x axes ranges
    if plotLimitsX:
        assert len(plotLimitsX) == len(axes1)
        for axis1,limits in zip(axes1,plotLimitsX):
            axis1.set_xlim(*limits)
    else:
        for axis1 in axes1:
            axis1.set_xlim(auto=True)
        print("Warning: No plot ranges given for X!")

    # The y axes ranges
    if plotLimitsY:
        limits=plotLimitsY
        assert len(plotLimitsY) == 2 
        assert len(plotLimitsY) == 2
        for axis1,axis2 in zip(axes1,axes2):
            axis1.set_ylim(*limits[0])
            axis2.set_ylim(*limits[1])

def _setAxesLabel(axes1,axes2,labelsXAxis,labelsYAxis,figure):
    """
    Sets the lables of the x and y axes using figure
    """
    # Labels for axes
    if labelsYAxis:
        assert len(axes1)==len(axes2)
        for axis1,axis2 in zip(axes1,axes2):
            axis1.set_ylabel(labelsYAxis[0])
            axis2.set_ylabel(labelsYAxis[1])

    if labelsXAxis:
        figure.text(0.45,0.005, labelsXAxis[0])

def _generateDiagonals(axes1,axes2,diagonalSize,distanceBetweenPlots,figure):
    """
    Generates the diagonals on x axis using a list of axes1 and 2 with
    the diagonalSize and the distance between subplots of the figure
    """
    # Generate Diagonals
    # arguments to pass to plot, just so we don't keep repeating them
    d = diagonalSize  # how big to make the diagonal lines in axes coordinates
    for index in range(0,len(axes1)) :
        if axes1[index] == axes1[-1]:
            break

        ax1=axes1[index]
        ax2=axes2[index+1]
        
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((1-d, 1+d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (1-d, 1+d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the right axes
        ax2.plot((-d,+d), (1-d,1+d), **kwargs)
        ax2.plot((-d,+d), (-d,+d), **kwargs)

    # Modify distance between plots
    figure.subplots_adjust(distanceBetweenPlots)

def _modifyVisability(axes1,axes2):
    """
    Sets the visibility of the axis for interrupted plot given a list of
    y-axis 1 and 2
    """
    # Hide spines and labels
    for axis1,axis2 in zip(axes1,axes2):
        axis1.yaxis.set_visible(False)
        removeSpines(axis1)
        axis2.yaxis.set_visible(False)
        removeSpines(axis2)

    # Show left y axis and right y axis
    axes1[0].yaxis.set_visible(True)
    axes1[0].spines['left'].set_visible(True)
    axes2[-1].yaxis.set_visible(True)
    axes2[-1].spines['right'].set_visible(True)

def _setXAxisTicks(axes1,majorTicksX,minorTicksX):
    """
    Sets the axis ticks for X given the list of axes to a FixedLocator
    majorTicksX and minorTicksX (as an array of floats)
    """
    if majorTicksX and minorTicksX:
        #for axis,majors,minors in zip(axes1,majorTicksX,minorTicksX):
        for axis in axes1:
            axis.xaxis.set_major_locator(ticker.FixedLocator(majorTicksX))
            axis.xaxis.set_minor_locator(ticker.FixedLocator(minorTicksX))
    

def setupPlot(numberOfXAxisInterruptions,plt,plotLimitsX=None,
        plotLimitsY=None,labelsYAxis="",labelsXAxis="",
        diagonalSize=.010, distanceBetweenPlots=0.03,
        majorTicksX=None,minorTicksX=None,
        majorTicksY=None,minorTicksY=None,
        tight=True):
    """
    This function does the setup for a Plot with interrupted XAxis. The number
    of interruptions can be specified via the *numberOfXAxisInterruptions* option


        :param plt: the pyplot object
        :param plotLimitsX: A list of tuples for each axis interruption
        :param plotLimitsY: A list of tuples for each axis interruption
        :param labelsYAxis: A list of length 2 with YAxis Labels
        :param labelsXAxis: Label for the x Axis. Will be placed approximatly in the center.
        :param digonalSize: Size of the diagonals at the axis interruption
        :param distanceBetweenPlots: distance between the Subplots
        :param majorTicksX: A list of lists each containing the fixed ticks location
        :param minorTicksX: A list of lists each containing the fixed ticks location
        :param majorTicksY: A list of lists each containing the fixed ticks location
        :param minorTicksY: A list of lists each containing the fixed ticks location

    Returns the figure and a list of length 2 for each y Axis
    
    """

    figure,axes1,axes2=_createFigureAndAxes(numberOfXAxisInterruptions,tight)
    _setAxesRanges(axes1,axes2,plotLimitsX,plotLimitsY)
    _setAxesLabel(axes1,axes2,labelsXAxis,labelsYAxis,figure)
    _generateDiagonals(axes1,axes2,diagonalSize,distanceBetweenPlots,figure)
    _modifyVisability(axes1,axes2)
    _setXAxisTicks(axes1,majorTicksX,minorTicksX)

    return figure,[axes1,axes2]

    
def plot(filenames,axes,labels,colors,
        yshift,figure,invert=False,
        peakWidth=5.0,lineWidthPlot=0.7,assignments=None,offsetAssignment=0.01):
    """
    Plot the spectra
        :param filenames: The filennames of the data
        :param axes: The axes to use in the plot (one for each subplot. See :ref setupPlot)
        :param labels: The lables for each plot
        :param colors: The colors for each plot
        :param yshift: Shift the data in y direction
        :param figure: The figure to be generated (See :ref setupPlot)
        :param invert: Select of data should be inverted
        :param peakWidth: For the assignments: Select the peak within ±peakWidth
        :param lineWidthPlot: The line width of the plots
        :param assignments: The assignments for the data (peaks)
        :param offsetAssignment: The vertical offset in axis coordinates

    Returns the plots and the assignmentPoints in axis coordinate
    """
    # Loop through all arguments for plotting
    plts=[]
    assignmentPoints=[]
    for i,filename in enumerate(filenames):
        # Import Data
        xdata,ydata=importData(filename)
        # Shift the data and invert it
        if invert:
            ydata*=-1
        ydata+=yshift
        # Now plot the data on each subplot
        plt=None
        for axis in axes:
            plt=axis.plot(xdata,ydata,linewidth=lineWidthPlot,
                    color=colors[i],
                    label=labels[i])

        # Compute the y value for each x value in the assignment 
        # These are in axis coordinates
        if assignments and all(assignments):
            for assignment in assignments[i]:
                # Assignment not none or nan
                if assignment:
                    peakWidth=peakWidth
                    yvalue,lowest=getPeakYForX(xdata,ydata,assignment,peakWidth)

                    # If the peak cannot be found
                    if not yvalue:
                        assignmentPoints.append((None,None,None))
                        continue

                    # Set the offset
                    if lowest:
                        offset=offsetAssignment*-1
                    else:
                        offset=offsetAssignment

                    for index,axis in zip(range(0,len(axes)),axes):
                        # Convert to axis coordinates
                        xAxisPoint,yAxisPoint=axis.transLimits.transform(
                                (assignment,yvalue))

                        # Only add point if within the axis view
                        if xAxisPoint >= 1.0 or xAxisPoint <= 0.0:
                            continue

                        else:
                            assignmentPoints.append(
                                    (index,xAxisPoint,yAxisPoint+offset))
                            break

                    # If no axis can be found
                    else:
                        assignmentPoints.append((None,None,None))

                # Assignment is None or NaN
                else:
                    assignmentPoints.append((None,None,None))

        plts+=plt
    return plts,assignmentPoints

def drawAssignments(pointsExp,pointsComp,axes,linewidth=0.3):
    for pointExp,pointComp in zip(pointsExp,pointsComp):
        axesIndexExp,xExp,yExp=pointExp
        axesIndexComp,xComp,yComp=pointComp
        if xExp and yExp:
            axis=axes[axesIndexExp]
            axis.annotate("", (xExp,yExp),
                    xytext=(xComp,yComp),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-",ls='--',
                        lw=linewidth))

