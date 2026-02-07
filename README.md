# MagTools

Tools for modeling magnetic source bodies


## Description

MagTools is a set of tools for modeling magnetized source bodies in Matlab. The software is organized in an object-oriented framework with several interacting classes.

Typical usage involves defining a magnetized object and a region over which to display the resulting field and then using functions like showBfieldContours or showBfieldVectors to visualize the resulting magnetic field. 

See examples.m for a detailed explanation and working examples. For even more examples, see reproduce_paper_figures.m. 


## Getting Started

MagTools depends on a library called BaseTools, which is typically located in a folder alongside the MagTools folder or under a subfolder called libs.

At the beginning of each session, the user must run setup.m to bring both MagTools and BaseTools (another library containing low-level functions) into the Matlab workspace.

The file called examples.m provides a basic introduction to using MagTools and includes lots of internal documentation.


## Disclaimer

This code is provided as-is, has been tested only very informally, and may not always behave as intended. It is actively under development and future versions may not be backward compatible. The authors do not guarantee accuracy or robustness.

Maintenance note: This repository is shared in the interest of Open Science. While you are free to use and adapt the code under the MIT License, we do not provide technical support, bug fixes, or guarantee future compatibility.
