## Accumulation of distances in `double`.

I think that this should be fine for sketch sizes up to 1e12 or so, because
this is the point at which any likely signal (say 0.001) is represented exactly
in a double.

However, for safety, we sum over each chunk independently, and then sum over
the chunk.

I think this is a good balance between speed and the saftey of something like
[fsum](https://code.activestate.com/recipes/393090/) or [Kahan
summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm)


