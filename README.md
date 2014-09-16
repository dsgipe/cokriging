cokriging
=======

A class used to create cokriging surrogate models
depencies: class Arr in array project

uses:
Cokriging is a statistical method used to combine sparse data sets found through expensive means, with more prevalent data determined from a source of lower confidence. For example, combining 100 inviscid CFD points with 10 viscous CFD points to build a model as if 100 viscous CFD points were run. These methods have been shown to significantly decrease cost for developing accurate meta models. As far as the author knows, there is no free or non propriatary  available cokriging software that can handle multiple dimensions. Currently this class can only do one dimension. But furthering this to N-dimensions is the final goal.

More information can be found by reading Forrester's papers on kriging and cokriging

To do: Replace pointers with more effective memory management class Arr, which is in the array project.
    Update to work with 2D data sets
