# Data Keys

Data keys are a hyphen delimited sequence of alpha numeric characters which describe a particular configuration of the billiard. An example of a `Data Key` might be something like `s1-b0-p2-t1--`. The breakdown for this arbitrary key is given below:

    - `b0` corresponds to the 0th `boundary parameter` (think the radius of a circular scatterer or location of polygon vertex). Of course, there may be multiple boundary parameters involved 
     
    - `p2` is a `partition` index which says that there are multiple instances (here 2 or 3 depending on whether you count the 0th partition) of this billiard configuration. Useful for partitioning an ensemble into smaller pieces.

    - `t1` could refer to the `total time` being a variable quantity. This is used with finite domain systems where particles that don't escape within a fixed time are not considered (one then wants to compare quantities which depend on the total time allowed). 

    - `s1` is a generic identifier which points to a particular results set which may be overlapping with other `systems`. This can happen when there are algorithmic changes to the execution of the billiard but all other parameters are the same. For example, we can simulate the same billiard `b0` over the same time interval `t0` with the same number of particles, with the only difference being that `s1` uses a rejection sampling algorithm while `s0` does not. 
     Because each `Study` is limited to one executable file, the occurrence of different `systems` will need to be handled from within  


From the above we see that what counts as a "configuration" is, in principle an endless list of combinations of all the different variables involved. Most likely is that we don't want to do anything too complicated but one never knows how expansive his exploration will be so it is better to leave the door open.


# Directives 

A `Directive` is an association between the name of a data file and the corresponding methods which can be used organise that data. Different data files depend on `Data Keys` in different ways. For instance, consider the following file list:

    1. `b0-p0--data.dat`
    2. `b0-p1--data.dat`
    3. `b0-p2--data.dat`
    4. `b1-p0--data.dat`
    5. `b1-p1--data.dat`

Suppose `data.dat` has dimensions `M x N` (though `1 x N` or `M x 1` can also be used. Inconsistent dimensions may also be possible). Then there are different `directives` or ways we combine, manage or transform the data. For instance, we might like to concatenate all unique `bx`'s over the different partitions:

    `cat b0-p*--data.dat > b0--data.dat`
    `cat b1-p*--data.dat > b1--data.dat`

Perhaps `data.dat` is actually a matrix of averages. Then we would like to perform a transformation before collapsing the different partitions into one file:

    `Rscript MatrixMeans.R *data.dat` 
    `cat b0-p*--data.dat > b0--data.dat`
    `cat b1-p*--data.dat > b1--data.dat`

However when possible we average data internally so as to save on steps in post processing. In fact, the vast majority of data can be handled by a simple collapsing of partitions. 

SUGGESTION: Is there anything other than collapsing partitions?




