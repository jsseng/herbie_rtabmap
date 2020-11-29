# Cached RTABMAP

This is a fork of RTABMAP that caches datafiles for fast initialization on startup.  The code assumes that you have built a map using RTABMAP already and would like to use the map just for localization.

### Notes

* The code currently supports only ORB features (256-bit binary descriptors)
* Benchmarking
  * Without the cache, startup time were about 10 seconds (on an Intel Xeon E5-1620).  With the cache, startup time are about 0.7 seconds.

### Usage

* Once you have built a map, build the cache files using:

```
rtabmap  --Mem/IncrementalMemory false --Mem/InitWMWithAllNodes true --Rtabmap/BuildCache true database.db
```

When this command is run, is will load the database file and build the cached data files.  You should see the new files: `vwdictionary.dat` and `flann.dat` created in your directory.  You can close RTABMAP once the database has loaded.

* To start RTABMAP with the cache files:

```
rtabmap  --Mem/IncrementalMemory false --Mem/InitWMWithAllNodes true --Rtabmap/UseCache true database.db
```


rtabmap ![Analytics](https://ga-beacon-279122.nn.r.appspot.com/UA-56986679-3/github-main?pixel) 
=======

[![RTAB-Map Logo](https://raw.githubusercontent.com/introlab/rtabmap/master/guilib/src/images/RTAB-Map100.png)](http://introlab.github.io/rtabmap)

[![Release][release-image]][releases]
[![License][license-image]][license]
Linux: [![Build Status](https://travis-ci.org/introlab/rtabmap.svg?branch=master)](https://travis-ci.org/introlab/rtabmap) Windows: [![Build status](https://ci.appveyor.com/api/projects/status/hr73xspix9oqa26h/branch/master?svg=true)](https://ci.appveyor.com/project/matlabbe/rtabmap/branch/master)

[release-image]: https://img.shields.io/badge/release-0.20.2-green.svg?style=flat
[releases]: https://github.com/introlab/rtabmap/releases

[license-image]: https://img.shields.io/badge/license-BSD-green.svg?style=flat
[license]: https://github.com/introlab/rtabmap/blob/master/LICENSE

RTAB-Map library and standalone application.

For more information, visit the [RTAB-Map's home page](http://introlab.github.io/rtabmap) or the [RTAB-Map's wiki](https://github.com/introlab/rtabmap/wiki).

To use RTAB-Map under ROS, visit the [rtabmap](http://wiki.ros.org/rtabmap) page on the ROS wiki.

### Acknowledgements
This project is supported by [IntRoLab - Intelligent / Interactive / Integrated / Interdisciplinary Robot Lab](https://introlab.3it.usherbrooke.ca/), Sherbrooke, Québec, Canada.

<a href="https://introlab.3it.usherbrooke.ca/">
<img src="https://github.com/introlab/16SoundsUSB/blob/master/images/IntRoLab.png" alt="IntRoLab" height="100">
</a>
