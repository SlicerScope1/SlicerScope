import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging


import SimpleITK as sitk
import sys
import numpy

import math


import os.path
from os.path import expanduser
homeDir = expanduser("~")


try:
    import openslide
except ImportError:
    slicer.util.pip_install("openslide-python")

try:
    import skimage.transform
except ImportError:
    slicer.util.pip_install("scilearn-image")


try:
    import h5py
except ImportError:
    slicer.util.pip_install("h5py")

try:
    import keras.models
except ImportError:
    slicer.util.pip_install("keras")



os.environ["CUDA_VISIBLE_DEVICES"] = "0"


#
# BigViewerModule
#

class BigViewerModule(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "BigViewerModule" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Big Viewer"]
    self.parent.dependencies = []
    self.parent.contributors = ["Yi Gao"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    BigViewer for viewing large image whose whole content is not able to be loaded into memory.
    """
    self.parent.acknowledgementText = """
    This file was developed by Yi Gao.
""" # replace with organization, grant and thanks.

#
# BigViewerModuleWidget
#

class BigViewerModuleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent = None):
    ScriptedLoadableModuleWidget.__init__(self, parent) # must have this. otherwise get error when loading this module


    print("+++++++++++++++++++++++++++++  BigViewerModuleWidget", flush = True)

    self.kerasModelSegmentNuclei = keras.models.load_model("/home/gaoyi/usr/work/project/forSlicerScope/SlicerScope/Modules/ServerSide/unet_nucleus_20200120_at20200126.hdf5")

    self.kerasModelSegmentColonGland = keras.models.load_model("/home/gaoyi/usr/work/project/forSlicerScope/SlicerScope/Modules/ServerSide/unet_glas_20200713.hdf5")



    # set to Red view only
    lm = slicer.app.layoutManager()
    lm.setLayout(6)

    #--------------------------------------------------------------------------------
    # Viewer variables
    self.topLeftX = 0
    self.topLeftY = 0

    self.patchSizeX = 0
    self.patchSizeY = 0

    self.MPP = 0.25 # um/pixel
    self.SliceThickness = 0.005 #mm

    self.BigRGBAImagePathname = ""

    self.BigRGBAImageNumberOfLevels = 0
    self.BigRGBAImageLevelToLoad = 0

    self.leftMouseButtonPos = (0, 0)

    self.WSISizesXAtAllLevels = []
    self.WSISizesYAtAllLevels = []

    self.H5FilePathname = ""
    self.h5FileLoaded = False

    print("I'm in 11111111111111111111111111111111111111111")


    # Interator
    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    interactor = redView.interactorStyle().GetInteractor()
    interactor.AddObserver(vtk.vtkCommand.RightButtonPressEvent, self.onRightButtonPressed)

    interactor.AddObserver(vtk.vtkCommand.LeftButtonPressEvent, self.onLeftButtonPressed)
    interactor.AddObserver(vtk.vtkCommand.LeftButtonReleaseEvent, self.onLeftButtonReleased)

    interactor.AddObserver(vtk.vtkCommand.MouseWheelForwardEvent, self.onMouseWheelForwardEvent)
    interactor.AddObserver(vtk.vtkCommand.MouseWheelBackwardEvent, self.onMouseWheelBackwardEvent)

  #
  # Customized mouse right button pressed event
  #
  def onRightButtonPressed(self, obj, event=None):
    print ('onRightButtonPressed............................', event)

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    style = redView.interactorStyle()
    interactor = style.GetInteractor()
    print(interactor.GetEventPosition())
    # Do something here



  def onLeftButtonPressed(self, obj, event=None):
    print ('onLeftButtonPressed............................', event)

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    style = redView.interactorStyle()
    interactor = style.GetInteractor()

    self.leftMouseButtonPos = interactor.GetEventPosition()
    print(self.leftMouseButtonPos, type(self.leftMouseButtonPos), type(self.leftMouseButtonPos[0]))
    # Do something here


  def onLeftButtonReleased(self, obj, event=None):
    print ('onLeftButtonReleased............................', event)

    layoutManager = slicer.app.layoutManager()
    redWidget = layoutManager.sliceWidget('Red')
    redView = redWidget.sliceView()
    style = redView.interactorStyle()
    interactor = style.GetInteractor()
    print(interactor.GetEventPosition())

    newPos = interactor.GetEventPosition()

    self.topLeftXSliderWidget.setValue(self.topLeftXSliderWidget.value + 10.0*(self.leftMouseButtonPos[0] - newPos[0]))
    self.topLeftYSliderWidget.setValue(self.topLeftYSliderWidget.value + 10.0*(newPos[1] - self.leftMouseButtonPos[1])) # x and y are different so the orders are different

    #self.topLeftYSliderWidget.update()



  def onMouseWheelForwardEvent(self, obj, event=None):
    print ('onMouseWheelForwardEvent............................')
    if not self.ObjectiveMagnificationSlicerWidget.isMaximized():
      #self.ObjectiveMagnificationSlicerWidget.setValue(self.ObjectiveMagnificationSlicerWidget.value + 2.0*self.ObjectiveMagnificationSlicerWidget.singleStep)

      # set the increment to be 1/2 of the current values accelerates
      # the zooming well when the mag is large. It also make the
      # zooming stable when the value is small. This is better than
      # setting the increament to a fixed value
      self.ObjectiveMagnificationSlicerWidget.setValue(self.ObjectiveMagnificationSlicerWidget.value + 0.5*self.ObjectiveMagnificationSlicerWidget.value)

      # If send this signal every time the slider changes, will trigger loading patch. This will make the zooming very slow
      #self.ObjectiveMagnificationSlicerWidget.valueChanged(self.ObjectiveMagnificationSlicerWidget.value)


  def onMouseWheelBackwardEvent(self, obj, event=None):
    print ('onMouseWheelBackwardEvent............................')
    if not self.ObjectiveMagnificationSlicerWidget.isMinimized():
      #self.ObjectiveMagnificationSlicerWidget.setValue(self.ObjectiveMagnificationSlicerWidget.value - 2.0*self.ObjectiveMagnificationSlicerWidget.singleStep)

      # set the increment to be 1/2 of the current values accelerates
      # the zooming well when the mag is large. It also make the
      # zooming stable when the value is small. This is better than
      # setting the increament to a fixed value
      self.ObjectiveMagnificationSlicerWidget.setValue(self.ObjectiveMagnificationSlicerWidget.value - 0.5*self.ObjectiveMagnificationSlicerWidget.value)

      # If send this signal every time the slider changes, will trigger loading patch. This will make the zooming very slow
      #self.ObjectiveMagnificationSlicerWidget.valueChanged(self.ObjectiveMagnificationSlicerWidget.value)

  def setup(self):

    print("I'm in 2222222222222222222222222222222")

    # In this function, we instantiate and connect widgets ...
    ScriptedLoadableModuleWidget.setup(self)

    self.logic = BigViewerModuleLogic()

    #--------------------------------------------------------------------------------
    # Parameters Area
    WSIParametersCollapsibleButton = ctk.ctkCollapsibleButton()
    WSIParametersCollapsibleButton.text = "BigViewer"
    self.layout.addWidget(WSIParametersCollapsibleButton)

    # Layout within the dummy collapsible button
    wsiParametersFormLayout = qt.QFormLayout(WSIParametersCollapsibleButton)

    #--------------------------------------------------------------------------------
    # BigRGBAImage filename. Note this only pick the file pathname, does not
    # load the volume.
    self.BigRGBAImageFileNameEditor = ctk.ctkPathLineEdit()
    self.BigRGBAImageFileNameEditor.setCurrentPath("")
    wsiParametersFormLayout.addRow("Select WSI: ", self.BigRGBAImageFileNameEditor)

    #--------------------------------------------------------------------------------
    # Load BigRGBAImage meta information to update UI
    self.loadWSIMetaInfoButton = qt.QPushButton("Load WSI")
    self.loadWSIMetaInfoButton.toolTip = "Load information from WSI to populate the module."
    self.loadWSIMetaInfoButton.enabled = True
    wsiParametersFormLayout.addRow(self.loadWSIMetaInfoButton)

    # #--------------------------------------------------------------------------------
    # # Load BigRGBAImage meta information. The meta info will help adjust the sliders max value etc.
    # self.loadBigRGBAImageButton = qt.QPushButton("Load WSI")
    # self.loadBigRGBAImageButton.toolTip = "Load the WSI"
    # self.loadBigRGBAImageButton.enabled = False
    # wsiParametersFormLayout.addRow(self.loadBigRGBAImageButton)

    #--------------------------------------------------------------------------------
    # Top left position in the level 0 of the WSI
    self.topLeftXSliderWidget = ctk.ctkSliderWidget()
    self.topLeftXSliderWidget.tracking = True

    self.topLeftXSliderWidget.singleStep = 10
    self.topLeftXSliderWidget.minimum = 0
    self.topLeftXSliderWidget.maximum = 100
    self.topLeftXSliderWidget.decimals = 0
    self.topLeftXSliderWidget.value = 0
    self.topLeftXSliderWidget.setToolTip("Top Left Corner, X-position.")
    self.topLeftXSliderWidget.enabled = False
    wsiParametersFormLayout.addRow("Top Left X", self.topLeftXSliderWidget)

    self.topLeftYSliderWidget = ctk.ctkSliderWidget()
    self.topLeftYSliderWidget.tracking = True
    self.topLeftYSliderWidget.singleStep = 10
    self.topLeftYSliderWidget.minimum = 0

    self.topLeftYSliderWidget.maximum = 100
    self.topLeftYSliderWidget.decimals = 0
    self.topLeftYSliderWidget.value = 0
    self.topLeftYSliderWidget.setToolTip("Top Left Corner, Y-position.")
    self.topLeftYSliderWidget.enabled = False
    wsiParametersFormLayout.addRow("Top Left Y", self.topLeftYSliderWidget)


    #--------------------------------------------------------------------------------
    # Level in the WSI
    self.ObjectiveMagnificationSlicerWidget = ctk.ctkSliderWidget()
    self.ObjectiveMagnificationSlicerWidget.tracking = True

    self.ObjectiveMagnificationSlicerWidget.singleStep = 0.1
    self.ObjectiveMagnificationSlicerWidget.minimum = 0
    self.ObjectiveMagnificationSlicerWidget.maximum = 100
    self.ObjectiveMagnificationSlicerWidget.decimals = 1
    self.ObjectiveMagnificationSlicerWidget.value = 1
    self.ObjectiveMagnificationSlicerWidget.setToolTip("Zooming")
    self.ObjectiveMagnificationSlicerWidget.enabled = False
    wsiParametersFormLayout.addRow("Zoom", self.ObjectiveMagnificationSlicerWidget)

    self.wsiLevelToLoad = 0


    #--------------------------------------------------------------------------------
    # Parameters Area
    H5ParametersCollapsibleButton = ctk.ctkCollapsibleButton()
    H5ParametersCollapsibleButton.text = "H5Parameters"
    self.layout.addWidget(H5ParametersCollapsibleButton)

    # Layout within the dummy collapsible button
    H5ParametersFormLayout = qt.QFormLayout(H5ParametersCollapsibleButton)

    #--------------------------------------------------------------------------------
    # Binary segmentation BigTiff file. Note this only pick the file pathname, does not
    # load the volume.
    self.H5FileFileNameEditor = ctk.ctkPathLineEdit()
    self.H5FileFileNameEditor.setCurrentPath("")
    H5ParametersFormLayout.addRow("Select H5 File: ", self.H5FileFileNameEditor)


    #--------------------------------------------------------------------------------
    # Load Label Image meta information. The meta info will help adjust the sliders max value etc.
    self.loadH5FileButton = qt.QPushButton("Load H5 File")
    self.loadH5FileButton.toolTip = "Load the H5 Image"
    self.loadH5FileButton.enabled = True
    H5ParametersFormLayout.addRow(self.loadH5FileButton)

    self.h5DatasetOption = qt.QComboBox()
    #self.h5DatasetOption.addItems(("asdf", "123123"))
    self.h5DatasetOption.enabled = False
    H5ParametersFormLayout.addRow(self.h5DatasetOption)

    #--------------------------------------------------------------------------------
    # check box to trigger taking screen shots for later use in tutorials
    self.extractHematoxylinOnFlyCheckBox = qt.QCheckBox()
    self.extractHematoxylinOnFlyCheckBox.checked = False
    self.extractHematoxylinOnFlyCheckBox.setToolTip("If checked, will extract Hematoxylin channel on the fly.")
    H5ParametersFormLayout.addRow("Realtime Extract Hematoxylin?", self.extractHematoxylinOnFlyCheckBox)


    #--------------------------------------------------------------------------------
    # Process Area
    ProcessPatchCollapsibleButton = ctk.ctkCollapsibleButton()
    ProcessPatchCollapsibleButton.text = "ProcessPatch"
    self.layout.addWidget(ProcessPatchCollapsibleButton)

    # Layout within the dummy collapsible button
    ProcessPatchFormLayout = qt.QFormLayout(ProcessPatchCollapsibleButton)

    # Button for segmenting nuclei
    self.decomposeStainButton = qt.QPushButton("Decompose Staining")
    self.decomposeStainButton.toolTip = "Decompose staining in this patch"
    self.decomposeStainButton.enabled = True
    ProcessPatchFormLayout.addRow(self.decomposeStainButton)

    # Button for segmenting nuclei
    self.segmentNucleiButton = qt.QPushButton("Segment Nuclei")
    self.segmentNucleiButton.toolTip = "Segment nuclei in this patch"
    self.segmentNucleiButton.enabled = True
    ProcessPatchFormLayout.addRow(self.segmentNucleiButton)

    # Button for gland detection
    self.detectGlandButton = qt.QPushButton("Segment Colon Gland")
    self.detectGlandButton.toolTip = "Detect gland figures in this patch"
    self.detectGlandButton.enabled = True
    ProcessPatchFormLayout.addRow(self.detectGlandButton)

    # # Button for mitosis detection
    # self.detectMitosisButton = qt.QPushButton("Detect Mitotic Figure")
    # self.detectMitosisButton.toolTip = "Detect mitotic figures in this patch"
    # self.detectMitosisButton.enabled = True
    # ProcessPatchFormLayout.addRow(self.detectMitosisButton)

    #--------------------------------------------------------------------------------
    # connections
    self.loadWSIMetaInfoButton.connect('clicked(bool)', self.onLoadWSIMetaInfoButton)

    #self.loadBigRGBAImageButton.connect('clicked(bool)', self.onLoadBigRGBAImageButton)
    self.loadH5FileButton.connect('clicked(bool)', self.onLoadH5FileButton)

    self.topLeftXSliderWidget.connect("valueChanged(double)", self.loadPatchFromBigRGBAImage)
    self.topLeftYSliderWidget.connect("valueChanged(double)", self.loadPatchFromBigRGBAImage)

    self.ObjectiveMagnificationSlicerWidget.connect("valueChanged(double)", self.onWSILevelChanged)

    self.segmentNucleiButton.connect('clicked(bool)', self.onSegmentNucleiButton)
    self.decomposeStainButton.connect('clicked(bool)', self.onDecomposeStainButton)
    self.detectGlandButton.connect('clicked(bool)', self.onDetectGlandButton)
    # connections
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






    print("I'm in 3333333333333333333333333333333333333333333333333333")






    # Add vertical spacer
    self.layout.addStretch(1)

  def setupForScalarPatch(self):
    imageSize = [self.patchSizeX, self.patchSizeY, 1]
    imageSpacing = [self.MPP/1000., self.MPP/1000., self.SliceThickness]
    voxelType = vtk.VTK_UNSIGNED_CHAR

    #--------------------------------------------------------------------------------
    # Allocate space and create node to store the RGB patch
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 1)

    # Create volume node
    IJKToRASDirectionMatrix = vtk.vtkMatrix4x4()
    IJKToRASDirectionMatrix.SetElement(0, 0, -1)
    IJKToRASDirectionMatrix.SetElement(1, 1, -1)
    # now input patch is using LPS (which does not make sense for
    # 2D image....... but in CLI for color decomposition, the ITK
    # image is by default LPS. so the output image is not aligned
    # with the color image (ras)

    volumeNode = slicer.vtkMRMLScalarVolumeNode()
    volumeNode.SetName("currentPatchGrayChannel")
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
    volumeNode.SetAndObserveImageData(imageData)
    slicer.mrmlScene.AddNode(volumeNode)

    # Add volume to scene
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()

    slicer.mrmlScene.AddNode(displayNode)
    colorNode = slicer.util.getNode('Green')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    volumeNode.CreateDefaultStorageNode()

    # store the loaded patch to widget member
    self.patchGrayVolumeNode = volumeNode

  def setupVolumeNodeToStoreRGBPatch(self):
    imageSize = [self.patchSizeX, self.patchSizeY, 1]
    imageSpacing = [self.MPP/1000., self.MPP/1000., self.SliceThickness]
    voxelType = vtk.VTK_UNSIGNED_CHAR

    #--------------------------------------------------------------------------------
    # Allocate space and create node to store the RGB patch
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 3)

    # Create volume node
    IJKToRASDirectionMatrix = vtk.vtkMatrix4x4()
    IJKToRASDirectionMatrix.SetElement(0, 0, -1)
    IJKToRASDirectionMatrix.SetElement(1, 1, -1)
    # now input patch is using LPS (which does not make sense for
    # 2D image....... but in CLI for color decomposition, the ITK
    # image is by default LPS. so the output image is not aligned
    # with the color image (ras)

    volumeNode = slicer.vtkMRMLVectorVolumeNode()
    volumeNode.SetName("currentPatchFromBigRGBAImage")
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
    volumeNode.SetAndObserveImageData(imageData)
    slicer.mrmlScene.AddNode(volumeNode)

    # Add volume to scene
    displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)

    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    volumeNode.CreateDefaultStorageNode()

    # store the loaded patch to widget member
    self.patchVolumeNode = volumeNode

    #--------------------------------------------------------------------------------
    # Set the patch volume to BigRGBAImage loader so it knows where to store
    # patch data

    selectionNode = slicer.app.applicationLogic().GetSelectionNode()

    # This will select the RGB patch to the bacgkround of the Red view 
    selectionNode.SetReferenceActiveVolumeID(self.patchVolumeNode.GetID())

    # This will select the label mask to the label of the Red
    # view. Without this line, the label node will be in slicer but
    # not selected. You need to manually pick everytime
    #selectionNode.SetReferenceActiveLabelVolumeID(self.patchLabelVolumeNode.GetID())
    slicer.app.applicationLogic().PropagateVolumeSelection(0)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  def setupForLoadingPatchFromH5File(self):

    if self.H5FileFileNameEditor.currentPath and os.path.isfile(self.H5FileFileNameEditor.currentPath):
      self.H5FilePathname = self.H5FileFileNameEditor.currentPath

      #--------------------------------------------------------------------------------
      # Read some meta info from H5File to populate slicer UI

      # self.H5FileLoader.SetH5FileFileName(self.H5FilePathname)
      # self.H5FileLoader.Initialization()

      # self.H5FileNumberOfLevels = self.H5FileLoader.GetH5FileLevel()
      # self.H5FileSizeX0 = self.H5FileLoader.GetSizeX0()
      # self.H5FileSizeY0 = self.H5FileLoader.GetSizeY0()
      # self.MPP = self.H5FileLoader.GetMPP()

      self.hdf5File = h5py.File(self.H5FilePathname, 'r')

      # #--------------------------------------------------------------------------------
      # # Read UI info to determine the patch size to be extracted
      # lm = slicer.app.layoutManager()
      # redWidget = lm.sliceWidget('Red')
      # redView = redWidget.sliceView()

      # self.patchSizeX = redView.width
      # self.patchSizeY = redView.height

      imageSize = [self.patchSizeX, self.patchSizeY, 1]
      imageSpacing = [self.MPP/1000., self.MPP/1000., self.SliceThickness]
      voxelType = vtk.VTK_UNSIGNED_CHAR

      #--------------------------------------------------------------------------------
      # Allocate space and create node to store the label patch
      labelImageData = vtk.vtkImageData()
      labelImageData.SetDimensions(imageSize)
      labelImageData.AllocateScalars(voxelType, 1)

      IJKToRASDirectionMatrix = vtk.vtkMatrix4x4()
      IJKToRASDirectionMatrix.SetElement(0, 0, -1)
      IJKToRASDirectionMatrix.SetElement(1, 1, -1)
      # now input patch is using LPS (which does not make sense for
      # 2D image....... but in CLI for color decomposition, the ITK
      # image is by default LPS. so the output image is not aligned
      # with the color image (ras)

      labelNode = slicer.vtkMRMLLabelMapVolumeNode()
      labelNode.SetName("currentPatchFromH5File-label")
      labelNode.SetSpacing(imageSpacing)
      labelNode.SetIJKToRASDirectionMatrix(IJKToRASDirectionMatrix)
      labelNode.SetAndObserveImageData(labelImageData)
      slicer.mrmlScene.AddNode(labelNode)

      a = slicer.util.array("currentPatchFromH5File-label")
      a[:] = 0

      # Add volume to scene
      labelDisplayNode = slicer.vtkMRMLLabelMapVolumeDisplayNode()
      slicer.mrmlScene.AddNode(labelDisplayNode)
      #labelColorNode = slicer.util.getNode('Random')
      labelColorNode = slicer.util.getNode('GenericColors')
      labelDisplayNode.SetAndObserveColorNodeID(labelColorNode.GetID())
      labelNode.SetAndObserveDisplayNodeID(labelDisplayNode.GetID())
      labelNode.CreateDefaultStorageNode()

      red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
      red_cn = red_logic.GetSliceCompositeNode()
      red_cn.SetLabelVolumeID(labelNode.GetID())

      # store the loaded label patch to widget member
      self.patchLabelVolumeNode = labelNode


      #--------------------------------------------------------------------------------
      # Set the patch volume to H5File loader so it knows where to store
      # patch data
      #self.H5FileLoader.SetOutputPatchImage(self.patchLabelVolumeNode.GetImageData())

      selectionNode = slicer.app.applicationLogic().GetSelectionNode()

      # This will select the RGB patch to the bacgkround of the Red view 
      #selectionNode.SetReferenceActiveVolumeID(self.patchVolumeNode.GetID())

      # This will select the label mask to the label of the Red
      # view. Without this line, the label node will be in slicer but
      # not selected. You need to manually pick everytime
      selectionNode.SetReferenceActiveLabelVolumeID(self.patchLabelVolumeNode.GetID())
      slicer.app.applicationLogic().PropagateVolumeSelection(0)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def cleanup(self):
    pass

  # def onSelect(self):
  #   self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()


  def onWSILevelChanged(self):
    ds = (self.objectiveMagnificationMax)/(self.ObjectiveMagnificationSlicerWidget.value)

    slide = openslide.OpenSlide(self.BigRGBAImagePathname)
    self.BigRGBAImageLevelToLoad = slide.get_best_level_for_downsample(ds)
    slide.close()

    #print(self.BigRGBAImageLevelToLoad)

    self.updateUIWidget()

    self.loadPatchFromBigRGBAImage()

  def enableAndInitUIWidget(self):
    #self.loadBigRGBAImageButton.enabled = True

    #--------------------------------------------------------------------------------
    # Set the max of slider for top-left corner
    # self.topLeftXSliderWidget.value = self.topLeftX
    # self.topLeftYSliderWidget.value = self.topLeftY

    self.topLeftXSliderWidget.value = 0
    self.topLeftYSliderWidget.value = 0

    self.topLeftXSliderWidget.maximum = self.WSISizesXAtAllLevels[0] - self.patchSizeX*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1
    self.topLeftYSliderWidget.maximum = self.WSISizesYAtAllLevels[0] - self.patchSizeY*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1

    self.topLeftXSliderWidget.enabled = True
    self.topLeftYSliderWidget.enabled = True

    #--------------------------------------------------------------------------------
    # Set the max of number of levels
    self.ObjectiveMagnificationSlicerWidget.maximum = self.objectiveMagnificationMax
    self.ObjectiveMagnificationSlicerWidget.minimum = self.objectiveMagnificationMin
    self.ObjectiveMagnificationSlicerWidget.enabled = True


  def updateUIWidget(self):
    # self.topLeftXSliderWidget.value = self.topLeftX
    # self.topLeftYSliderWidget.value = self.topLeftY

    self.topLeftXSliderWidget.maximum = self.WSISizesXAtAllLevels[0] - self.patchSizeX*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1
    self.topLeftYSliderWidget.maximum = self.WSISizesYAtAllLevels[0] - self.patchSizeY*self.level_downsamples[self.BigRGBAImageLevelToLoad] - 1

    #self.ObjectiveMagnificationSlicerWidget.maximum = self.BigRGBAImageNumberOfLevels - 1


  def getRegionFromFileAsRGBNumpyArray(self, WSIName, level, topX0, topY0, width, height):
    # I currently do not have range check. May add later
    slide = openslide.OpenSlide(WSIName)

    thisTilePilIm = slide.read_region((topX0, topY0), level, (width, height))
    if thisTilePilIm.mode != "RGB":
        thisTilePilIm = thisTilePilIm.convert("RGB")

    imRGBNdarray = numpy.asarray(thisTilePilIm)

    slide.close()

    return imRGBNdarray

  def loadPatchFromBigRGBAImage(self):
    #--------------------------------------------------------------------------------
    # Load RGB image
    #
    # The image to view is of the size (self.patchSizeX,
    # self.patchSizeY) on the screen. However, i need to compute how
    # at the current zooming level (self.BigRGBAImageLevelToLoad), how
    # big a region i need to load, and then zoom to (self.patchSizeX,
    # self.patchSizeY)
    #
    #

    #print(self.BigRGBAImageLevelToLoad)
    # print(int(self.topLeftXSliderWidget.value))
    # print(int(self.topLeftYSliderWidget.value))
    # print(int(self.patchSizeX))
    # print(int(self.patchSizeY))


    magAtThisLevel = float(self.slideInfo["objectiveMagnification"])/float(self.level_downsamples[self.BigRGBAImageLevelToLoad])

    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    # print(self.ObjectiveMagnificationSlicerWidget.value)
    ratio = magAtThisLevel/float(self.ObjectiveMagnificationSlicerWidget.value)
    # print(ratio)
    # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    sizeToLoadX = round(self.patchSizeX*ratio)
    sizeToLoadY = round(self.patchSizeY*ratio)

    self.topLeftX = int(self.topLeftXSliderWidget.value)

    #print("RGB requested size = ", (sizeToLoadY, sizeToLoadX))
    q = self.getRegionFromFileAsRGBNumpyArray(self.BigRGBAImagePathname, self.BigRGBAImageLevelToLoad, int(self.topLeftXSliderWidget.value), int(self.topLeftYSliderWidget.value), int(sizeToLoadX), int(sizeToLoadY))

    # print(q.max())

    p = skimage.transform.resize(q, (int(self.patchSizeY), int(self.patchSizeX)), preserve_range=True, mode='wrap', order=0, anti_aliasing=False)
    #print("RGB Got size = ", q.shape, "resize to ", p.shape, q.shape[0]/p.shape[0], q.shape[1]/p.shape[1])
    # print(p.max())

    #a = slicer.util.array('currentPatchFromBigRGBAImage')
    volumeNode = slicer.util.getNode('currentPatchFromBigRGBAImage')
    a = slicer.util.arrayFromVolume(volumeNode)

    # resizing is needed for python (not for CLI version) coz we have
    # to fill the numpy array here
    a[:] = p

    resolutionNow = self.MPP/1000.0*float(self.slideInfo["objectiveMagnification"])/float(self.ObjectiveMagnificationSlicerWidget.value)
    imageSpacing = [resolutionNow, resolutionNow, self.SliceThickness]

    n = slicer.util.getNode('currentPatchFromBigRGBAImage')
    n.SetSpacing(imageSpacing)

    n.GetImageData().Modified()

    #--------------------------------------------------------------------------------
    # Load label image
    if self.h5FileLoaded:
      group = self.hdf5File[self.h5DatasetOption.currentText]

      dsetName = "data" + str(self.BigRGBAImageLevelToLoad)
      dset = group[dsetName]

      start0 = int(round(self.topLeftYSliderWidget.value/self.level_downsamples[self.BigRGBAImageLevelToLoad]))
      start1 = int(round(self.topLeftXSliderWidget.value/self.level_downsamples[self.BigRGBAImageLevelToLoad]))
      size0 = int(round(self.patchSizeY*ratio))
      size1 = int(round(self.patchSizeX*ratio))

      # print(dsetName)
      # print((start0, start1))
      #print(size1)

      #print("requested size = ", (size0, size1))

      labelArrayInLevel = dset[start0:(start0 + size0), start1:(start1 + size1)]
      #print("max label labelArrayInLevel = " + str(labelArrayInLevel.max()))

      #---------------------------------------------------------------
      # When at the borader, when requested region is to the right or
      # below the read image, the openslide returns an array the same
      # size as the requested size. But h5 only return the intersected
      # region. So the return image size is smaller than the requested
      # size. If not padded first, the subsequent resize will stretch
      # the label image
      if labelArrayInLevel.shape != (size0, size1):
        tmp = numpy.zeros((size0, size1))
        tmp[:labelArrayInLevel.shape[0], :labelArrayInLevel.shape[1]] = labelArrayInLevel
        labelArrayInLevel = tmp
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      labelArray = skimage.transform.resize(labelArrayInLevel, (int(self.patchSizeY), int(self.patchSizeX)), preserve_range=True, mode='wrap', order=0, anti_aliasing=False)
      #print("got size = ", labelArrayInLevel.shape, "resize to = ", labelArray.shape, "ratio = ", labelArrayInLevel.shape[0]/labelArray.shape[0], labelArrayInLevel.shape[1]/labelArray.shape[1])
      #print("max label labelArray = " + str(labelArray.max()))

      #labelArray = dset[int(self.topLeftYSliderWidget.value):int(self.topLeftYSliderWidget.value + self.patchSizeY), int(self.topLeftXSliderWidget.value):int(self.topLeftXSliderWidget.value + self.patchSizeX)]

      # print(self.topLeftXSliderWidget.value, self.topLeftXSliderWidget.value, self.patchSizeX, self.patchSizeY)
      # print(labelArray.max())

      # labelArray = numpy.zeros((self.patchSizeY, self.patchSizeX), dtype='uint8')
      # labelArray[int(self.patchSizeY/3):int(2*self.patchSizeY/3), int(self.patchSizeX/3):int(2*self.patchSizeX/3)] = 1

      n = slicer.util.getNode('currentPatchFromH5File-label')
      a = slicer.util.array('currentPatchFromH5File-label')
      a[:] = labelArray

      n.SetSpacing(imageSpacing)
      n.GetImageData().Modified()

      b = slicer.util.array('currentPatchFromH5File-label')

      self.patchLabelVolumeNode.GetImageData().Modified()

    #--------------------------------------------------------------------------------
    # If need
    if self.extractHematoxylinOnFlyCheckBox.checked:
      parameters = {}
      parameters['inputVolume'] = self.patchVolumeNode.GetID()
      parameters['outputVolume'] = self.patchGrayVolumeNode.GetID()
      slicer.cli.run( slicer.modules.colordecomposition, None, parameters, wait_for_completion=True )



    #--------------------------------------------------------------------------------
    # magnitude = vtk.vtkImageMagnitude()
    # magnitude.SetInputData(self.patchVolumeNode.GetImageData())
    # magnitude.Update()

    # grayNode = slicer.vtkMRMLScalarVolumeNode()
    # grayNode.SetImageDataConnection(magnitude.GetOutputPort())
    # grayNode.SetName("tempGrayNode")

    # labelNode = slicer.vtkMRMLLabelMapVolumeNode()
    # logic = BigViewerModuleLogic()
    # logic.run(grayNode, labelNode, 100)

    # slicer.mrmlScene.AddNode(labelNode)

    # # Add volume to scene
    # displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()

    # slicer.mrmlScene.AddNode(displayNode)
    # colorNode = slicer.util.getNode('FreeSurferLabels')
    # displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    # labelNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    # labelNode.CreateDefaultStorageNode()

    # red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
    # red_cn = red_logic.GetSliceCompositeNode()
    # red_cn.SetLabelVolumeID(labelNode.GetID())


    # magnitude = vtk.vtkImageMagnitude()
    # magnitude.SetInputData(self.patchVolumeNode.GetImageData())
    # magnitude.Update()

    # grayNode = slicer.vtkMRMLScalarVolumeNode()
    # grayNode.SetImageDataConnection(magnitude.GetOutputPort())
    # grayNode.SetName("tempGrayNode")

    # slicer.mrmlScene.AddNode(grayNode)

    # displayNode = slicer.vtkMRMLVectorVolumeDisplayNode()

    # slicer.mrmlScene.AddNode(displayNode)
    # colorNode = slicer.util.getNode('FreeSurferLabels')
    # displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    # grayNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    # grayNode.CreateDefaultStorageNode()

    # red_logic = slicer.app.layoutManager().sliceWidget("Red").sliceLogic()
    # red_cn = red_logic.GetSliceCompositeNode()
    # red_cn.SetForegroundVolumeID(grayNode.GetID())


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # #bgrdNode.SetIJKToRASDirectionMatrix(fMat)
    # slicer.mrmlScene.AddNode(bgrdNode)
    # bgrdVolID = bgrdNode.GetID()  
    # red_cn.SetForegroundVolumeID(fgrdVolID)
    # red_cn.SetBackgroundVolumeID(bgrdVolID)
    # red_cn.SetForegroundOpacity(1)   



    # bgrdNode.SetImageDataConnection(magnitude.GetOutputPort())
    # bgrdNode.SetName(bgrdName)
    # #bgrdNode.SetIJKToRASDirectionMatrix(fMat)
    # slicer.mrmlScene.AddNode(bgrdNode)
    # bgrdVolID = bgrdNode.GetID()  
    # red_cn.SetForegroundVolumeID(fgrdVolID)
    # red_cn.SetBackgroundVolumeID(bgrdVolID)
    # red_cn.SetForegroundOpacity(1)   


#    enableScreenshotsFlag = self.extractHematoxylinOnFlyCheckBox.checked



    # redViewWidget = slicer.app.layoutManager().sliceWidget("Red")
    # redView = redViewWidget.sliceView()
    # redView.

    lm = slicer.app.layoutManager()
    redWidget = lm.sliceWidget('Red')
    redView = redWidget.sliceView()

    redController = redWidget.sliceController()
    redController.fitSliceToBackground()
    #redView.resetFocalPoint()
    redView.forceRender()

    #    r = slicer.app.layoutManager().sliceWidget("Red").sliceController()
#    r.fitSliceToBackground()


    # wti = vtk.vtkWindowToImageFilter()
    # wti.SetInput(redView.renderWindow())
    # wti.Update()
    # v = vtk.vtkImageViewer()
    # v.SetColorWindow(255)
    # v.SetColorLevel(128)
    # v.SetInputConnection(wti.GetOutputPort())
    # v.Render()

    #slicer.app.processEvents()
    #qt.QApplication.processEvents()

  def determinPatchSizeByViewerSize(self):
    #--------------------------------------------------------------------------------
    # Read UI info to determine the patch size to be extracted
    lm = slicer.app.layoutManager()
    redWidget = lm.sliceWidget('Red')
    redView = redWidget.sliceView()

    self.patchSizeX = int(redView.width/2)
    self.patchSizeY = int(redView.height/2)

    return

  def getMetaInfoFromWSI(self):
    logic = BigViewerModuleLogic()
    #self.slideInfo = hptk.getSlideInfo(self.BigRGBAImagePathname)
    self.slideInfo = logic.getSlideInfo(self.BigRGBAImagePathname)

    self.BigRGBAImageNumberOfLevels = self.slideInfo["level_count"]
    self.level_downsamples = self.slideInfo["level_downsamples"]
    self.MPP = self.slideInfo["mpp"][0]

    for it in range(self.BigRGBAImageNumberOfLevels):
      self.WSISizesXAtAllLevels.append(self.slideInfo["level_dimensions"][it][0])
      self.WSISizesYAtAllLevels.append(self.slideInfo["level_dimensions"][it][1])

    self.objectiveMagnificationMax = self.slideInfo["objectiveMagnification"]
    self.objectiveMagnificationMin = self.slideInfo["objectiveMagnification"]/self.level_downsamples[-1]/3.0 # so to include more

    # This will trigger the onWSILevelChanged in which the
    # loadPatchFromBigRGBAImage will be called, where the
    # currentPatchGrayChannel node will be filled--- but it's not
    # allocated yet. So this should not be set here.
    #
    # self.ObjectiveMagnificationSlicerWidget.value = self.objectiveMagnificationMax 

    # print(self.objectiveMagnificationMax)
    # print(self.objectiveMagnificationMin)

    return

  def onLoadWSIMetaInfoButton(self):
    self.determinPatchSizeByViewerSize()

    #--------------------------------------------------------------------------------
    # Read meta infrom from WSI image
    #if self.BigRGBAImageFileNameEditor.currentPath and os.path.isfile(self.BigRGBAImageFileNameEditor.currentPath):
    if not os.path.isfile(self.BigRGBAImageFileNameEditor.currentPath):
      print("File not exist")
      return

    self.BigRGBAImagePathname = self.BigRGBAImageFileNameEditor.currentPath
    self.getMetaInfoFromWSI()

    self.setupVolumeNodeToStoreRGBPatch()
    self.setupForScalarPatch()

    # Setting the magnification slide bar will automatically trigger
    # the onWSILevelChanged in which will call the
    # loadPatchFromBigRGBAImage
    #self.ObjectiveMagnificationSlicerWidget.value = self.objectiveMagnificationMax
    self.ObjectiveMagnificationSlicerWidget.value = self.objectiveMagnificationMin

    # Setting the magnification slide bar will automatically trigger
    # the onWSILevelChanged in which will call the
    # loadPatchFromBigRGBAImage. So no need to do again
    # self.loadPatchFromBigRGBAImage()

    self.enableAndInitUIWidget()

  # def onLoadBigRGBAImageButton(self):
  #   self.setupVolumeNodeToStoreRGBPatch()
  #   self.loadPatchFromBigRGBAImage()

  #   self.enableAndInitUIWidget()



  def onDecomposeStainButton(self):
    parameters = {}
    parameters['inputVolume'] = self.patchVolumeNode.GetID()
    parameters['outputVolume'] = self.patchGrayVolumeNode.GetID()
    slicer.cli.run( slicer.modules.colordecomposition, None, parameters, wait_for_completion=True )

    volumeNode = slicer.util.getNode('currentPatchGrayChannel') #self.patchGrayVolumeNode
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage') #self.patchVolumeNode
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)


  def onDetectMitosisButton(self):
    pass

  def onDetectGlandButton(self):
    self.logic.processSegmentGland(self.patchVolumeNode, self.patchGrayVolumeNode, self.kerasModelSegmentColonGland)

    volumeNode = slicer.util.getNode('currentPatchGrayChannel')
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage')
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)

    return

  def onSegmentNucleiButton(self):
    self.logic.processSegmentNuclei(self.patchVolumeNode, self.patchGrayVolumeNode, self.kerasModelSegmentNuclei)

    volumeNode = slicer.util.getNode('currentPatchGrayChannel')
    slicer.util.setSliceViewerLayers(foreground=volumeNode)
    volumeNode1 = slicer.util.getNode('currentPatchFromBigRGBAImage')
    slicer.util.setSliceViewerLayers(background = volumeNode1)

    slicer.util.setSliceViewerLayers(foregroundOpacity=0.4)

  def onLoadH5FileButton(self):
    if not self.patchVolumeNode:
      return

    self.setupForLoadingPatchFromH5File()

    self.h5FileLoaded = True

    f = h5py.File(self.H5FilePathname, 'r')

    self.groupNamesInH5 = []

    for key in f.keys():
        if isinstance(f[key], h5py.Group):
          # node is a dataset
          print("Group ", key)
          self.groupNamesInH5.append(key)
          self.h5DatasetOption.addItem(key)
        elif isinstance(f[key], h5py.Dataset):
          #self.groupNamesInH5.append(key)
          print("Dataset ", key)

    f.close()

    self.h5DatasetOption.enabled = True



#
# BigViewerModuleLogic
#

class BigViewerModuleLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True


  def getSlideInfo(self, svsPathname):
    slide = openslide.OpenSlide(svsPathname)

    slideInfo = {'slidePathname':svsPathname,
                 'level_count':slide.level_count,
                 'level_dimensions': slide.level_dimensions,
                 'level_downsamples': slide.level_downsamples,
                 'mpp': [float(slide.properties['openslide.mpp-x']), float(slide.properties['openslide.mpp-y'])],
                 'objectiveMagnification': float(slide.properties['openslide.objective-power'])}
    slide.close()

    return slideInfo



  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, imageThreshold, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : threshold, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('BigViewerModuleTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True

  def segment2DRGBPatchBatch(self, model, inputImagePatchBatchArray):
    # This segments a list of RGB 2D patches. The input
    # inputImagePatchArray is a (#batches, patchSideLen, patchSideLen,
    # numChannel) numpy array. Next this needs to be reshaped to
    # (#batches, patchSideLen, patchSideLen, #channels)
    # where and #channels=1 in this case

    sz = inputImagePatchBatchArray.shape

    inputImagePatchBatchArray /= 255.0

    testBatchSize = sz[0]
    results = model.predict(inputImagePatchBatchArray, testBatchSize, verbose=1)

    outputSegBatchArray = results[:, :, :, :]

    return outputSegBatchArray


  def segment2DRGBImageRandomSampleDividePrior(self, model, imageArray, patchSideLen = 64, numPatchSampleFactor = 10, batch_size = 1, num_segmetnation_classes = 3):
#def segment2DRGBImageRandomSampleDividePrior(model, imageArray, patchSideLen = 64, numPatchSampleFactor = 10, batch_size = 1, num_segmetnation_classes = 3):
    sz = imageArray.shape
    numChannel = 3 # for RGB

    print("+++++++++++++++++++", sz)

    #assert(sz[0] >= patchSideLen and sz[1] >= patchSideLen and sz[2] == 3),"Image shape must be >= " + str(patchSideLen) + "-cubed."
    if sz[2] != numChannel:
        print("Only process RGB image")
        exit(-1)

    # the number of random patches is s.t. on average, each pixel is
    # sampled numPatchSampleFactor times. Default is 10
    numPatchSample = math.ceil((sz[0]/patchSideLen)*(sz[1]/patchSideLen)*numPatchSampleFactor)


    print("numPatchSample = ", numPatchSample)

    # this saves the segmentation result
    segArray = numpy.zeros((sz[0], sz[1], num_segmetnation_classes), dtype=numpy.float32)
    priorImage = numpy.zeros((sz[0], sz[1]), dtype=numpy.float32)

    patchShape = (patchSideLen, patchSideLen, numChannel)
    imagePatchBatch = numpy.zeros((batch_size, patchShape[0], patchShape[1], numChannel), dtype=numpy.float32)

    for itPatch in range(0, numPatchSample, batch_size):

        allPatchTopLeftX = numpy.random.randint(0, sz[0] - patchShape[0], size = batch_size)
        allPatchTopLeftY = numpy.random.randint(0, sz[1] - patchShape[1], size = batch_size)

        for itBatch in range(batch_size):
            thisTopLeftX = allPatchTopLeftX[itBatch]
            thisTopLeftY = allPatchTopLeftY[itBatch]

            imagePatchBatch[itBatch, :, :, :] = imageArray[thisTopLeftX:(thisTopLeftX + patchShape[0]), thisTopLeftY:(thisTopLeftY + patchShape[1]), :]

        segBatch = self.segment2DRGBPatchBatch(model, imagePatchBatch)

        for itBatch in range(batch_size):
            thisTopLeftX = allPatchTopLeftX[itBatch]
            thisTopLeftY = allPatchTopLeftY[itBatch]

            segArray[thisTopLeftX:(thisTopLeftX + patchShape[0]), thisTopLeftY:(thisTopLeftY + patchShape[1]), :] += segBatch[itBatch, :, :, :]
            priorImage[thisTopLeftX:(thisTopLeftX + patchShape[0]), thisTopLeftY:(thisTopLeftY + patchShape[1])] += numpy.ones((patchShape[0], patchShape[1]))

    for it in range(num_segmetnation_classes):
        segArray[:, :, it] /= (priorImage + numpy.finfo(numpy.float32).eps)
        segArray[:, :, it] *= 100

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # segArray contains multiple channels of output. The 1-st is the
    # corresponding output for the 1st object.
    outputSegArrayOfObject1 = segArray[:, :, 1]

    return outputSegArrayOfObject1

  def processSegmentGland(self, inputVolume, outputVolume, kerasModel):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param outputVolume: thresholding result
    :param imageThreshold: values above/below this threshold will be set to 0
    :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
    :param showResult: show output volume in slice viewers
    """

    if not inputVolume or not outputVolume:
      raise ValueError("Input or output volume is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

    ################################################################################
    # keras seg
    imgArray = slicer.util.array("currentPatchFromBigRGBAImage")
    imgArray = imgArray[0, :, :, :]

    testImgSeg = self.segment2DRGBImageRandomSampleDividePrior(model = kerasModel, imageArray = imgArray.astype(numpy.float), patchSideLen = 64, numPatchSampleFactor = 5, batch_size = 10)
    print(testImgSeg.max(), flush = True)

    testImgSeg = testImgSeg[numpy.newaxis, :, :]

    # volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    # outputVolume.CreateDefaultDisplayNodes()
    slicer.util.updateVolumeFromArray(outputVolume, testImgSeg)
    #setSliceViewerLayers(background=volumeNode)
    slicer.util.setSliceViewerLayers(foreground=outputVolume)


    # keras seg, end
    ################################################################################

    stopTime = time.time()
    logging.info('Processing completed in {0:.2f} seconds'.format(stopTime-startTime))



  def processSegmentNuclei(self, inputVolume, outputVolume, kerasModel):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputVolume: volume to be thresholded
    :param outputVolume: thresholding result
    :param imageThreshold: values above/below this threshold will be set to 0
    :param invert: if True then values above the threshold will be set to 0, otherwise values below are set to 0
    :param showResult: show output volume in slice viewers
    """

    if not inputVolume or not outputVolume:
      raise ValueError("Input or output volume is invalid")

    import time
    startTime = time.time()
    logging.info('Processing started')

    # # Compute the thresholded output volume using the "Threshold Scalar Volume" CLI module
    # cliParams = {
    #   'InputVolume': inputVolume.GetID(),
    #   'OutputVolume': outputVolume.GetID(),
    #   'ThresholdValue' : imageThreshold,
    #   'ThresholdType' : 'Above' if invert else 'Below'
    #   }
    # cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True, update_display=showResult)
    # # We don't need the CLI module node anymore, remove it to not clutter the scene with it
    # slicer.mrmlScene.RemoveNode(cliNode)


    ################################################################################
    # keras seg
    imgArray = slicer.util.array("currentPatchFromBigRGBAImage")
    imgArray = imgArray[0, :, :, :]

    testImgSeg = self.segment2DRGBImageRandomSampleDividePrior(model = kerasModel, imageArray = imgArray.astype(numpy.float), patchSideLen = 64, numPatchSampleFactor = 5, batch_size = 10)
    print(testImgSeg.max(), flush = True)

    testImgSeg = testImgSeg[numpy.newaxis, :, :]

    # volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
    # outputVolume.CreateDefaultDisplayNodes()
    slicer.util.updateVolumeFromArray(outputVolume, testImgSeg)
    #setSliceViewerLayers(background=volumeNode)
    slicer.util.setSliceViewerLayers(foreground=outputVolume)


    # keras seg, end
    ################################################################################

    stopTime = time.time()
    logging.info('Processing completed in {0:.2f} seconds'.format(stopTime-startTime))


class BigViewerModuleTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_BigViewerModule1()

  def test_BigViewerModule1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = BigViewerModuleLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
