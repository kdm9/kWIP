# `kWIP` 2.0

# Dependencies

System dependencies are:

- ``liblz4``
- ``zlib``
- ``libhdf5``
- ``gsl``

I.e.:

```bash
    sudo aptitude install libhdf5-dev zlib1g-dev liblz4-dev libgsl-dev
```

In addition, these are bundled:

- ``c-blosc`` and the C ``hdf5`` filter
- ``xxhash``
- ``greatest``
- ``kseq.h``
