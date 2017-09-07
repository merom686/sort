#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <functional>
#include <iterator>
#include <cmath>
#include <climits>

template <class I, class C = std::less<>>
void insertion_sort(I first, I last, C c = C()) {
    for (I it = first; ++it < last;) {
        auto v0 = *(it - 1);
        const auto v1 = *it;
        if (c(v1, v0)) {
            I i = it;
            do {
                *i-- = v0;
            } while (i > first && c(v1, v0 = *(i - 1)));
            *i = v1;
        }
    }
}

template <class I, class C = std::less<>>
void shell_sort(I first, I last, C c = C()) {
    size_t n = last - first;
    size_t h = 1;
    while (h * 9 * 3 < n) h = h * 9 + 1;

    for (; h > 0; h /= 9) {
        for (size_t k = 0; k < h; k++) {
            for (I it = first + k; (it += h) < last;) {
                auto v0 = *(it - h);
                const auto v1 = *it;
                if (c(v1, v0)) {
                    I i = it;
                    do {
                        *i = v0;
                        i -= h;
                    } while (i - h >= first && c(v1, v0 = *(i - h)));
                    *i = v1;
                }
            }
        }
    }
}

template <class I, class C = std::less<>>
void bubble_sort(I first, I last, C c = C()) {
    for (;;) {
        bool b = false;
        for (I it = first + 1; it < last; ++it) {
            if (c(*it, *(it - 1))) std::swap(*it, *(it - 1)), b = true;
        }
        if (!b) break;
    }
}

template <class I, class C = std::less<>>
void comb_sort(I first, I last, C c = C()) {
    for (size_t h = last - first;;) {
        if (h > 1) h = h * 10 / 13;
        bool b = false;
        for (I it = first + h; it < last; ++it) {
            if (c(*it, *(it - h))) std::swap(*it, *(it - h)), b = true;
        }
        if (h <= 1 && !b) break;
    }
}

#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("no-tree-pre")
#endif
template <class I, class C>
size_t larger_child(I p, C c, size_t j) {
    size_t k = j * 2;
    if (c(p[k], p[k + 1])) k++;
    return k;
}

template <class I, class C>
void down_heap(I p, C c, size_t i, size_t j) {
    auto t = p[i]; p[i] = p[j];
    size_t m = i / 2;
    while (j < m) {
        size_t k = larger_child(p, c, j);
        if (!c(t, p[k])) break;
        p[j] = p[k];
        j = k;
    }
    p[j] = t;
}

template <class I, class C = std::less<>>
void heap_sort(I first, I last, C c = C()) {
    size_t n = last - first;
    if (n <= 48) {
        insertion_sort(first, last, c);
        return;
    }
    I p = first - 1;

    // nÇ™ãÙêîÇÃÇ∆Ç´ p[n]Ç™íPì∆ÇÃç≈ëÂÇ≈Ç†Ç¡ÇƒÇÕÇ»ÇÁÇ»Ç¢ p[n]ÇèúÇ¢Çƒç\íz
    // ç≈ëÂílÇéÊÇËèoÇµÇΩÇÁp[n]Ç∆åä∑Ç∑ÇÈÇÃÇ≈ÇªÇÃéûì_Ç≈ïÅí ÇÃìÒï™ÉqÅ[ÉvÇ…Ç»ÇÈ
    size_t m = (n - 1) / 2;
    if (n % 2 == 0) {
        if (c(p[n / 2], p[n])) std::swap(p[n / 2], p[n]);
    }
    for (size_t i = m; i != 0; i--) {
        size_t k = larger_child(p, c, i);
        auto t = p[i];
        if (!c(t, p[k])) continue;
        p[i] = p[k];

        size_t j = k;
        while (j <= m) {
            k = larger_child(p, c, j);
            if (!c(t, p[k])) break;
            p[j] = p[k];
            j = k;
        }
        p[j] = t;
    }
    if (n % 2 == 0) {
        down_heap(p, c, n, 1);
    }
    for (size_t i = n - (n % 2 == 0) - 1; i > 2; i -= 2) {
        down_heap(p, c, i + 0, 2 + c(p[2], p[3]));
        down_heap(p, c, i + 1, 1);
    }
    auto t = p[3]; p[3] = p[1];
    if (c(p[2], t)) p[1] = p[2], p[2] = t; else p[1] = t;
}
#ifdef __GNUC__
#pragma GCC pop_options
#endif

template <class I, class P>
I partition(I first, I last, P p) {
    for (I i0 = first, i1 = last;;) {
        while (i0 < i1 && p(*i0)) i0++;
        do {
            if (i0 == i1) return i0;
        } while (!p(*--i1));
        std::swap(*i0++, *i1);
    }
}

template <class I, class C = std::less<>>
void quick_sort(I first, I last, C c = C()) {
    size_t n = last - first;
    if (n <= 64) {
        insertion_sort(first, last, c);
        return;
    }
    auto p = *(first + n / 2);
    I i0 = ::partition(first, last, [&](const auto& v) { return c(v, p); });
    quick_sort(first, i0, c);
    I i1 = ::partition(i0, last, [&](const auto& v) { return !c(p, v); });
    quick_sort(i1, last, c);

    //std::swap(*(first + n / 2), *(last - 1));
    //auto p = *(last - 1);
    //I it = ::partition(first, last - 1, [&](const auto& v) { return c(v, p); });
    //if (it < last - 1) std::swap(*it, *(last - 1));
    //quick_sort(first, it, c);
    //quick_sort(it + 1, last, c);
}

template <class I, class C>
void intro_sort_impl(I first, I last, int m, C c = C()) {
    size_t n = last - first;
    if (n <= 48) {
        insertion_sort(first, last, c);
        return;
    }
    if (m-- <= 0 || n < 1 << 12) {
        heap_sort(first, last, c);
        return;
    }
    std::swap(*(first + n / 2), *(last - 1));
    auto p = *(last - 1);
    I it = ::partition(first, last - 1, [&](const auto& v) { return c(v, p); });
    if (it < last - 1) std::swap(*it, *(last - 1));
    intro_sort_impl(first, it, m, c);
    intro_sort_impl(it + 1, last, m, c);

    //auto p = *(first + n / 2);
    //I i0 = ::partition(first, last, [&](const auto& v) { return c(v, p); });
    //if (first == i0) m = 0; else intro_sort_impl(first, i0, m, c);
    //intro_sort_impl(i0, last, m, c);
}

template <class I, class C = std::less<>>
void intro_sort(I first, I last, C c = C()) {
    intro_sort_impl(first, last, (int)std::log2(last - first), c);
}

template <class I, class T, class C>
void merge_sort_impl(I first, I last, T *temp, C c = C()) {
    size_t n = last - first;
    if (n <= 64) {
        insertion_sort(first, last, c);
        return;
    }
    size_t m = n / 2;
    I mid = first + m;
    merge_sort_impl(first, mid, temp, c);
    merge_sort_impl(mid, last, temp, c);

    std::copy(first, mid, temp);
    T *temp_last = temp + m, *i0 = temp;
    I it = first, i1 = mid;
    for (;;) {
        if (c(*i0, *i1)) {
            *it++ = *i0++;
            if (i0 == temp_last) break;
        } else {
            *it++ = *i1++;
            if (i1 == last) {
                std::copy(i0, temp_last, it);
                break;
            }
        }
    }
}

template <class I, class C = std::less<>>
void merge_sort(I first, I last, C c = C()) {
    using T = typename std::iterator_traits<I>::value_type;
    size_t n = last - first;
    T *temp = nullptr;
    if (n > 64) temp = (T *)std::malloc(n / 2 * sizeof(T));
    merge_sort_impl(first, last, temp, c);
    if (n > 64) free(temp);
}

template <class I, class T, class C = std::less<>>
inline I lower_bound(I i0, I i1, T t, C c = C()) {
    while (i0 < i1) {
        I i = i0 + (i1 - i0) / 2;
        if (!c(*i, t)) i1 = i; else i0 = i + 1;
    }
    return i0;
}

template <class I, class T, class C = std::less<>>
inline I upper_bound(I i0, I i1, T t, C c = C()) {
    while (i0 < i1) {
        I i = i0 + (i1 - i0) / 2;
        if (c(t, *i)) i1 = i; else i0 = i + 1;
    }
    return i0;
}

template <class I>
void reverse(I first, I last) {
    while (first < --last) std::swap(*first++, *last);
}

template <class I>
void rotate(I first, I mid, I last) {
    if (first == mid || mid == last) return;
    ::reverse(first, mid);
    ::reverse(mid, last);
    ::reverse(first, last);
}

template <class I, class C>
void merge_in_place(I first, I mid, I last, size_t n0, size_t n1, C c) {
    if (n0 == 0 || n1 == 0) return;
    if (n0 + n1 <= 40) {
        insertion_sort(first, last, c);
        return;
    }
    if (!c(*mid, *(mid - 1))) return;

    I i0, i1;
    size_t m0, m1;
    if (n0 >= n1) {
        // pà»â∫|pà»è„|pñ¢ñû|pà»è„
        m0 = n0 / 2;
        i0 = first + m0;
        i1 = ::lower_bound(mid, last, *i0, c);
        m1 = i1 - mid;
    } else {
        // pà»â∫|pÇÊÇËè„|pà»â∫|pà»è„
        m1 = (n1 - 1) / 2;
        i1 = mid + m1;
        i0 = ::upper_bound(first, mid, *i1, c);
        m0 = i0 - first;
    }
    // rotateÇ≈ìôÇµÇ¢óvëfÇí«Ç¢âzÇ∑Ç±Ç∆Ç™Ç»Ç¢ÇÊÇ§Ç…Ç∑ÇÈ
    ::rotate(i0, mid, i1);

    merge_in_place(first, i0, i0 + m1, m0, m1, c);
    merge_in_place(i0 + m1, i1, last, n0 - m0, n1 - m1, c);
}

template <class I, class C = std::less<>>
void in_place_merge_sort(I first, I last, C c = C()) {
    size_t n = last - first;
    if (n <= 64) {
        insertion_sort(first, last, c);
        return;
    }
    size_t m = n / 2;
    I mid = first + m;
    in_place_merge_sort(first, mid, c);
    in_place_merge_sort(mid, last, c);
    merge_in_place(first, mid, last, m, n - m, c);
}

template <class I, class C = std::less<>>
void radix_sort(I first, I last, C) {
    using T = typename std::iterator_traits<I>::value_type;
    static_assert(std::is_integral<T>::value, "");
    constexpr int L = 8, M = 1 << L, K = (sizeof(T) * CHAR_BIT + L - 1) / L;

    size_t n = last - first;
    T *p[2] = { &*first, (T *)std::malloc(n * sizeof(T)) };
    size_t m[M];

    for (int k = 0; k < K; k++) {
        std::fill_n(m, M, 0);
        for (size_t i = 0; i < n; i++) {
            m[(p[0][i] >> (L * k) & M - 1)]++;
        }
        for (int j = M - 1; j > 0; j--) {
            m[j] = m[j - 1];
        }
        m[0] = 0;
        for (int j = 2; j < M; j++) {
            m[j] += m[j - 1];
        }
        for (size_t i = 0; i < n; i++) {
            T t = p[0][i];
            p[1][m[t >> (L * k) & M - 1]++] = t;
        }
        std::swap(p[0], p[1]);
    }
    if (K % 2 == 1) {
        std::swap(p[0], p[1]);
        std::copy_n(p[1], n, p[0]);
    }
    if (std::is_signed<T>::value) {
        T *q = ::upper_bound(p[0], p[0] + n, 0, std::greater<>());
        ::rotate(p[0], q, p[0] + n);
    }
    std::free(p[1]);
}

template <class I, class C = std::less<>>
void std_heap_sort(I first, I last, C c = C()) {
    std::partial_sort(first, last, last, c);
}

void print(std::vector<int> v, size_t n) {
    for (size_t i = 0; i < n; i++) {
        std::cout << v[i] << ' ';
    }
    std::cout << '\n';
}

template <class F>
void bench(F f) {
    constexpr size_t N = 1000000, M = 500;
    static std::vector<int> v(N), u(N);

    using Duration = std::chrono::nanoseconds;
    static std::vector<Duration::rep> d(M);

    std::mt19937 rnd;
    for (size_t n = 1; n <= N; n *= 10) {
        size_t l = std::min(N * 5 / n, (size_t)M);
        for (size_t k = 0; k < l; k++) {
            for (size_t i = 0; i < n; i++) {
                v[i] = rnd(); //(rnd() >> (32 - 10)) - (int)i;
            }
            int64_t s = std::accumulate(v.begin(), v.begin() + n, 0LL);
            std::copy_n(v.begin(), n, u.begin());

            auto t0 = std::chrono::high_resolution_clock::now();
            f(v.begin(), v.begin() + n, std::less<>());
            auto t1 = std::chrono::high_resolution_clock::now();
            (t1 - t0).count();
            d[k] = std::chrono::duration_cast<Duration>(t1 - t0).count();

            if (!std::is_sorted(v.begin(), v.begin() + n) || s != std::accumulate(v.begin(), v.begin() + n, 0LL)) {
                print(v, n);
                std::sort(u.begin(), u.begin() + n);
                print(u, n);
                std::cout << "failed" << std::endl;
            }
        }
        std::sort(d.begin(), d.begin() + l);
        l = l * 2 / 3;
        std::cout << std::setw(9) << std::right << std::accumulate(d.begin(), d.begin() + l, (Duration::rep)0) / l << std::endl;
    }
}

int main() {
    using T = std::vector<int>::iterator;

    bench(heap_sort<T>);
    bench(intro_sort<T>);
    bench(quick_sort<T>);
    bench(merge_sort<T>);
    bench(in_place_merge_sort<T>);
    //bench(insertion_sort<T>);
    bench(shell_sort<T>);
    //bench(bubble_sort<T>);
    bench(comb_sort<T>);
    bench(radix_sort<T>);

    std::cout << "- STL -" << std::endl;
    bench(std::sort<T, std::less<>>);
    bench(std::stable_sort<T, std::less<>>);
    bench(std_heap_sort<T>);

    return 0;
}
