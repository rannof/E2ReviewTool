# E2ReviewTool
### ElarmS Review Tool - Review ElarmS E2 log files
Elarms is an Earthquake Early Warning System (EEWS), developed at Berkeley Seismological Lab (http://seismo.berkeley.edu/).
E2ReviewTool is aimed to help in reviewing ElarmS E2 module performenses, comparing results to a catalog file.
The tool allows a visual inspection of the events location and parameters including their evolution with time.
Created by Ran Novitsky Nof (ran.nof@gmail.com), 2015 @ BSL

#### NOTE: 
This is a beta version, and is still under development. Please report of any bug found.

#### DEPENDENCIES:
python (tested on 2.7) modules:

-   numpy - http://numpy.org
-   matplotlib - http://matplotlib.org
-   PyQt4 - http://www.riverbankcomputing.com
-   pyproj - https://pypi.python.org/pypi/pyproj
-   obspy - http://obspy.org

#### USAGE:
<pre>
E2ReviewTool.py [-h] [-i file [file ...]] [-r file [file ...]]
                       [-b bounds bounds bounds bounds]

optional arguments:
  -h, --help            show this help message and exit
  -i file [file ...]    input E2 log file(s)
  -r file [file ...]    input reference csv file csv - comma separated value
                        in GII DB format
                        [epiid,ml,typ,lat,lon,abs_ot_d,abs_ot_t,refDI]. epiid:
                        Event Catalog ID; ml: Local Magnitude; typ: Event type
                        code; lat: Latitude; lon: Longitude; abs_ot_d: Origin
                        time Date; abs_ot_t: Origin time time; refDI: E2
                        reference ID (or None)
  -b bounds bounds bounds bounds
                        Region bounding box (west east south north)
</pre>
#### LICENSE:
  E2ReviewTool is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


