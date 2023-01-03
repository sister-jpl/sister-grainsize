# SISTER Snow grain size PGE Documentation

## Description
The L2A snow grain size PGE takes as input surface reflectance and calculates snow grain size using the method of Nolin and Dozier (2000). Snow grain size is modeled as a function of scaled band area centered at the 1030 nm ice absorption feature:![](./examples/prisma_snow_spectrum.png)![](./examples/prisma_grainsize_example1.png)Nolin, A. W., & Dozier, J. (2000).A hyperspectral method for remotely sensing the grain size of snow.Remote sensing of Environment, 74(2), 207-216.[doi.org/10.1016/S0034-4257(00)00111-5](https://doi.org/10.1016/S0034-4257(00)00111-5)
## PGE Arguments

In addition to required MAAP job submission arguments the L2A snow grain size PGE also takes the following argument(s):


|Argument| Type |  Description | Default|
|---|---|---|---|
| l2a_granule| string |L2A corrected reflectance dataset granule URL| -|
| l2a_frac_cover| string |L2A fractional cover granule URL| -|
| frac_threshold| float |Fractional cover threshold| 0.6|
## Ouputs## Algorithm registration

## Exampple