## Accumulation of distances in `double`.

I think that this should be fine for sketch sizes up to 10G or so, because this
is the point at which any likely signal (say 0.001) is represented exactly in a
double. This was calculated as $2^53 / 255 * 0.001$, i.e. precision of double
divided by largest possible distance between `uint8_t`s, times the smallest
distance we are likely to care about.

However, for safety, we sum over each chunk independently, and then sum over
the chunk. I think this is a good balance between speed and the saftey of
something like [fsum](https://code.activestate.com/recipes/393090/) or [Kahan
summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm)

# Distance calculation is broken

We aren't calculating distances properly. We must fix this.

We will need to calculate the norms first (L1 or L2, or weighted L1/L2), then
element-wise calculate the normed values before calculating distances

Could also calulate L1 & L2 norms, and write into counts file, as part of
counting.
