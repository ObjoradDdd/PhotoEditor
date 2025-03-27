#define STB_IMAGE_IMPLEMENTATION
#include "libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "libs/stb_image_write.h"
#include <iostream>
#include <vector>
#include <cstdint>
#include <x86intrin.h>
#include <immintrin.h>
#include <cstdio>
#include <ctime>

using namespace std;


unsigned char posterize(unsigned char value, int levels = 5) {
    float step = 255.0f / (levels - 1);
    return (unsigned char) (int(int(value) / step) * step);
}

vector<unsigned char> toSepia(unsigned char R, unsigned char G, unsigned char B) {
    vector<unsigned char> rgb(3);
    rgb[0] = (unsigned char) min((int(R) * 0.393 + int(G) * 0.769 + int(B) * 0.189), 255.0);
    rgb[1] = (unsigned char) min((int(R) * 0.349 + int(G) * 0.686 + int(B) * 0.168), 255.0);
    rgb[2] = (unsigned char) min((int(R) * 0.272 + int(G) * 0.534 + int(B) * 0.131), 255.0);

    return rgb;
}

double sequentially(int repeat = 1, int code = 4) {
    int width, height, channels;
    double solarizationTime = 0, sepiaTime = 0, pasteurizationTime = 0;


    unsigned char *sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width, &height,
                                     &channels, 0);
    unsigned char *Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width,
                                            &height, &channels, 0);
    unsigned char *Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width,
                                              &height, &channels, 0);

    if (code == 1 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
            for (int i = 0; i < width * height * channels; ++i) {
                if (Solarization[i] > 128) {
                    Solarization[i] = (255 - Solarization[i]);
                }
            }
        }
        currentTime = clock() - currentTime;
        solarizationTime = (double) currentTime / CLOCKS_PER_SEC * 1000;
        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Solarization.png", width, height, channels,
                            Solarization, 0)) {
            stbi_image_free(Solarization);
        }
    }

    if (code == 2 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
            for (int i = 0; i < width * height * channels; i++) {
                Pasteurization[i] = posterize(Pasteurization[i]);
            }
        }
        currentTime = clock() - currentTime;
        pasteurizationTime = (double) currentTime / CLOCKS_PER_SEC * 1000;
        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Pasteurization.png", width, height, channels,
                            Pasteurization, 0)) {
            std::cout << "Pasteurization error" << std::endl;
            stbi_image_free(Pasteurization);
        }
    }

    if (code == 3 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
            for (int i = 0; i < width * height * channels - 3; i = i + 3) {
                vector<unsigned char> rgb = toSepia(sepia[i], sepia[i + 1], sepia[i + 2]);
                sepia[i] = rgb[0];
                sepia[i + 1] = rgb[1];
                sepia[i + 2] = rgb[2];
            }
        }
        currentTime = clock() - currentTime;
        sepiaTime = (double) currentTime / CLOCKS_PER_SEC * 1000;
        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Sepia.png", width, height, channels, sepia,
                            0)) {
            std::cout << "Sepia error" << std::endl;
            stbi_image_free(sepia);
        }
    }
    return sepiaTime + solarizationTime + pasteurizationTime;
}


double openMp(int repeat = 1, int code = 4) {
    int width, height, channels;
    double solarizationTime = 0, sepiaTime = 0, pasteurizationTime = 0;


    unsigned char *sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/openmp/Original.png", &width, &height,
                                     &channels, 0);
    unsigned char *Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/openmp/Original.png", &width,
                                            &height, &channels, 0);
    unsigned char *Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/openmp/Original.png", &width,
                                              &height, &channels, 0);

    if (code == 1 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
#pragma omp parallel for
            for (int i = 0; i < width * height * channels; ++i) {
                if (Solarization[i] > 128) {
                    Solarization[i] = 255 - Solarization[i];
                }
            }
        }
        currentTime = clock() - currentTime;
        solarizationTime = (double) currentTime / CLOCKS_PER_SEC * 1000;
        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/openmp/Solarization.png", width, height, channels,
                            Solarization, 0)) {
            stbi_image_free(Solarization);
        }
    }

    if (code == 2 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
#pragma omp parallel for
            for (int i = 0; i < width * height * channels; i++) {
                Pasteurization[i] = posterize(Pasteurization[i]);
            }
        }
        currentTime = clock() - currentTime;
        pasteurizationTime = (double) currentTime / CLOCKS_PER_SEC * 1000;
        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/openmp/Pasteurization.png", width, height, channels,
                            Pasteurization, 0)) {
            std::cout << "Pasteurization error" << std::endl;
            stbi_image_free(Pasteurization);
        }
    }

    if (code == 3 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
#pragma omp parallel for
            for (int i = 0; i < width * height * channels - 3; i = i + 3) {
                vector<unsigned char> rgb = toSepia(sepia[i], sepia[i + 1], sepia[i + 2]);
                sepia[i] = rgb[0];
                sepia[i + 1] = rgb[1];
                sepia[i + 2] = rgb[2];
            }
        }
        currentTime = clock() - currentTime;
        sepiaTime = (double) currentTime / CLOCKS_PER_SEC * 1000;
        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/openmp/Sepia.png", width, height, channels, sepia,
                            0)) {
            std::cout << "Sepia error" << std::endl;
            stbi_image_free(sepia);
        }
    }
    return sepiaTime + solarizationTime + pasteurizationTime;
}

double vectorize(int repeat = 1, int code = 4) {
    int width, height, channels;
    unsigned char *Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/vectorize/Original.png", &width, &height,
                                            &channels, 0);
    unsigned char *sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/vectorize/Original.png", &width, &height,
                                     &channels, 0);
    unsigned char *Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/vectorize/Original.png", &width, &height,
                                              &channels, 0);


    double solarizationTime = 0, sepiaTime = 0, pasteurizationTime = 0;
    if (code == 1 || code == 4) {
        clock_t currentTime = clock();
        for (int g = 0; g < repeat; ++g) {
            __m256i threshold16 = _mm256_set1_epi8(128);
            __m256i full255_16 = _mm256_set1_epi8(255);

            for (int i = 0; i < width * channels * height / 32 * 32 - 32; i += 32) {
                __m256i bytes = _mm256_load_si256((__m256i const *) &Solarization[i]);

                __m256i subtracted = _mm256_sub_epi8(full255_16, bytes);

                __m256i mask = _mm256_cmpeq_epi8(_mm256_max_epu8(bytes, threshold16), bytes);

                __m256i blended16 = _mm256_blendv_epi8(bytes, subtracted, mask);

                _mm256_storeu_si256((__m256i *) &Solarization[i], blended16);
            }
        }
        currentTime = clock() - currentTime;
        solarizationTime = (double) currentTime / CLOCKS_PER_SEC * 1000;

        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/vectorize/Solarization.png", width, height, channels,
                            Solarization, 0)) {
            stbi_image_free(Solarization);
        }
    }


    if (code == 2 || code == 4) {
        clock_t currentTime = clock();
        float step = 255.0f / 4.0f;
        __m256 step_vec = _mm256_set1_ps(step);
        for (int g = 0; g < repeat; ++g) {
            for (int i = 0; i + 8 < width * channels * height; i += 8) {
                __m256 pixels = _mm256_cvtepi32_ps(
                    _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i *) &Pasteurization[i]))
                );

                __m256 div = _mm256_div_ps(pixels, step_vec);
                __m256 floor = _mm256_floor_ps(div);

                __m256 result = _mm256_mul_ps(floor, step_vec);

                __m256i int32_result = _mm256_cvtps_epi32(result);

                __m128i low = _mm256_extractf128_si256(int32_result, 0);
                __m128i high = _mm256_extractf128_si256(int32_result, 1);

                __m128i packed16 = _mm_packs_epi32(low, high);

                __m128i packed8 = _mm_packus_epi16(packed16, packed16);

                _mm_storel_epi64((__m128i *) &Pasteurization[i], packed8);
            }
        }
        currentTime = clock() - currentTime;
        pasteurizationTime = (double) currentTime / CLOCKS_PER_SEC * 1000;

        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/vectorize/Pasteurization.png", width, height, channels,
                            Pasteurization, 0)) {
            stbi_image_free(Pasteurization);
        }
    }

    if (code == 3 || code == 4) {
        clock_t currentTime = clock();
        const int pixels = width * height;
        for (int g = 0; g < repeat; ++g) {
            const __m128 r_coef = _mm_setr_ps(0.393f, 0.349f, 0.272f, 0.0f);
            const __m128 g_coef = _mm_setr_ps(0.769f, 0.686f, 0.534f, 0.0f);
            const __m128 b_coef = _mm_setr_ps(0.189f, 0.168f, 0.131f, 0.0f);
            const __m128 max_val = _mm_set1_ps(255.0f);


            for (int i = 0; i < pixels; i++) {
                int idx = i * 3;

                unsigned char R = sepia[idx];
                unsigned char G = sepia[idx + 1];
                unsigned char B = sepia[idx + 2];

                __m128 r_vec = _mm_set1_ps((float) R);
                __m128 g_vec = _mm_set1_ps((float) G);
                __m128 b_vec = _mm_set1_ps((float) B);

                __m128 r_result = _mm_mul_ps(r_vec, r_coef);
                __m128 g_result = _mm_mul_ps(g_vec, g_coef);
                __m128 b_result = _mm_mul_ps(b_vec, b_coef);

                __m128 sum = _mm_add_ps(_mm_add_ps(r_result, g_result), b_result);

                __m128 clamped = _mm_min_ps(sum, max_val);

                float results[4];
                _mm_storeu_ps(results, clamped);

                sepia[idx] = (unsigned char) int(results[0]);
                sepia[idx + 1] = (unsigned char) int(results[1]);
                sepia[idx + 2] = (unsigned char) int(results[2]);
            }
        }
        currentTime = clock() - currentTime;
        sepiaTime = (double) currentTime / CLOCKS_PER_SEC * 1000;

        if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/vectorize/Sepia.png", width, height, channels,
                            sepia, 0)) {
            stbi_image_free(sepia);
        }
    }

    return pasteurizationTime + sepiaTime + solarizationTime;
}

int main() {
    printf("vectorize: %f milliseconds\n", vectorize(1, 3));
    printf("sequentially: %f milliseconds\n", sequentially(1, 3));
    printf("openMp: %f milliseconds\n", openMp(1, 4));
}
