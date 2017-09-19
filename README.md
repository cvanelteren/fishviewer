# Welcome to the fish page
Fast 3D time series visualization tool initially designed for the zebrafish dataset[Danio rerio].

This is a placeholder
Current features:

- 3D slicing
- Time series scrubbing
- Single cell time series selection
- ROI selection

# Dependencies
- numpy
- vispy
- matplotlib
- scipy
- sklearn

# Current Keybindings
Please note, most of the keybindings work on an active vispy canvas displaying
scatter data

- Control + shift + left click opens time data for the point
 - Due to a bug in OpenGL this only works when viewed from a flattened view
- Control + shift + right click adds time series to current active plot
- R, starts rotation
- Spacebar, will start autotime scroll
- Shift + mouse drag will shift center point of camera
- Scroll wheel will zoom
- Mouse drag will pan the camera
# Latest version showcase
![pageIntoGif][Video/Zebrafishviewer_Casper_van_Elteren.mp4]
# Old version showcase
![pageIntoGif](Video/pageIntro.gif)
