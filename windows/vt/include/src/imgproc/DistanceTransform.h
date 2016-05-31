#pragma once

namespace DistanceTransformInternal
{
    template <class T>
    void distanceTransformCore(T * const p, const int stride, const int n, int * const v, T * const z, T * const f)
    {
        for (int q = 0; q < n; q++)
        {
            f[q] = p[q*stride];
        }

        int k = 0;
        v[0] = 0;
        z[0] = ElTraits<T>::TypeMin();
        z[1] = ElTraits<T>::TypeMax();

        for (int q = 1; q < n; q++)
        {
    label:
            T s = ((f[q]+q*q) - (f[v[k]] + v[k]*v[k])) / (2*(q - v[k]));

            if (s <= z[k])
            {
                k--;
                goto label;
            }

            k++;
            v[k] = q;
            z[k] = s;
            z[k+1] = ElTraits<T>::TypeMax();
        }

        k = 0;

        for (int q = 0; q < n; q++)
        {
            while (z[k+1] < q)
            {
                k++;
            }

            p[q*stride] = (q-v[k])*(q-v[k]) + f[v[k]];
        }
    }
}

template<class T>
HRESULT distanceTransform(CTypedImg<T> & img)
{
    HRESULT hr = S_OK;

    const int w = img.Width(), h = img.Height();

    const int n = VtMax(w, h);
    vector<int> v;
    vector<T> z;
    vector<T> f;

    if ((hr = v.resize(n)) != S_OK)
        return hr;
    if ((hr = z.resize(n+1)) != S_OK)
        return hr;
    if ((hr = f.resize(n)) != S_OK)
        return hr;

    for (int y = 0; y < h; y++)
    {
        DistanceTransformInternal::distanceTransformCore(&img(0, y), 1, w, &v[0], &z[0], &f[0]);
    }

    for (int x = 0; x < w; x++)
    {
        DistanceTransformInternal::distanceTransformCore(&img(x, 0), img.StrideBytes()/sizeof(T), h, &v[0], &z[0], &f[0]);
    }

    return hr;
}

template<class T>
HRESULT distanceTransform3d(CTypedImg<T> & img)
{
    HRESULT hr = S_OK;

    const int w = img.Width(), h = img.Height(), d = img.Bands();

    const int n = VtMax(w, VtMax(h, d));
    vector<int> v;
    vector<T> z;
    vector<T> f;

    if ((hr = v.resize(n)) != S_OK)
        return hr;
    if ((hr = z.resize(n+1)) != S_OK)
        return hr;
    if ((hr = f.resize(n)) != S_OK)
        return hr;

    for (int pz = 0; pz < d; pz++)
    {
        for (int py = 0; py < h; py++)
        {
            DistanceTransformInternal::distanceTransformCore(&img(0, py, pz), d, w, &v[0], &z[0], &f[0]);
        }
    }

    for (int pz = 0; pz < d; pz++)
    {
        for (int px = 0; px < w; px++)
        {
            DistanceTransformInternal::distanceTransformCore(&img(px, 0, pz), w*d, h, &v[0], &z[0], &f[0]);
        }
    }

    for (int py = 0; py < h; py++)
    {
        for (int px = 0; px < w; px++)
        {
            DistanceTransformInternal::distanceTransformCore(&img(px, py, 0), 1, d, &v[0], &z[0], &f[0]);
        }
    }

    return hr;
}
