The purpose of PhyloBits is to be a lightweight Julia package that can (a) read/write
phylogenetic trees, without having to load a large list of dependencies and resulting 
slow computation times, and (b) provide additional useful basic functions used by 
BioGeoJulia/PhyBEARS, many of them replicating functionality originally available 
in the R package BioGeoBEARS (https://github.com/nmatzke/BioGeoBEARS). 

PhyloBits contains code from PhyloNetworks (https://github.com/crsl4/PhyloNetworks.jl) 
by Claudia Solis-Lemus and Cecile Ane, specifically the code for 
reading/writing phylogenetic trees, and the associated types.

auxiliary.jl
descriptive.jl
manipulateNet.jl
readwrite.jl
types.jl


The PhyloNetworks.jl package is licensed under the MIT "Expat" License:

> Copyright (c) 2014-2018: Claudia Solis-Lemus and Cecile Ane.
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
>
