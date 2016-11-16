# DVOS
Imperial College, Space Plasmas, Time Series Analysis library

This library provides data object handling for time series data
with emphasis on Space Physics data with standard metadata.

Data objects comprise sequences of multidimensional arrays of elements,
so a single scalar value would have a sequence length of 1 with 1 dimension of size 1.

Elements may be one of the following data types: 
double, integer, text (DvString), time (DvTime), Time Interval (DvEvent).
A data object will hold data of only one type.

Data objects may have lists of metadata attached, and metadata is also stored as a DVOS object.
Metadata may itself have metadata attached, resulting in a self similar hierarchy.

All common mathematical operations are supported on the data, and metadata checking
is performed to ensure operations are allowed (e.g. units, coordinate frames).
Objects resulting from mathematical operations will have relevant metadata derived
from the input objects.

Joining onto a common timeline is performed on the fly for analysis,
but explicit joining methods are also supported.
