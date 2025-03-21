#define STB_IMAGE_IMPLEMENTATION
#include "libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "libs/stb_image_write.h"
#include <iostream>
#include <vector>
#include <cstdint>
#include <x86intrin.h>
#include <immintrin.h>
#include <stdio.h>
#include <time.h>

using namespace std;

unsigned char solarize(unsigned char value) {
    if ((int) value > 128) {
        return (unsigned char) (255 - (int)value);
    }
    return value;
}

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

double sequentially() {
    int width, height, channels;

    unsigned char *sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width, &height,
                                     &channels, 0);
    unsigned char *Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width,
                                            &height, &channels, 0);
    unsigned char *Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width,
                                              &height, &channels, 0);

    clock_t start = clock();

    for (int i = 0; i < width * height * channels; i++) {
        Solarization[i] = solarize(Solarization[i]);
    }

    for (int i = 0; i < width * height * channels; i++) {
        Pasteurization[i] = posterize(Pasteurization[i]);
    }

    for (int i = 0; i < width * height * channels - 3; i = i + 3) {
        vector<unsigned char> rgb = toSepia(sepia[i], sepia[i + 1], sepia[i + 2]);
        sepia[i] = rgb[0];
        sepia[i + 1] = rgb[1];
        sepia[i + 2] = rgb[2];
    }

    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Sepia.png", width, height, channels, sepia, 0)) {
        std::cout << "Sepia error" << std::endl;
        stbi_image_free(sepia);
    }
    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Pasteurization.png", width, height, channels,
                        Pasteurization, 0)) {
        std::cout << "Pasteurization error" << std::endl;
        stbi_image_free(Pasteurization);
    }
    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Solarization.png", width, height, channels,
                        Solarization, 0)) {
        
        stbi_image_free(Solarization);
    }
    return seconds;
}

double openMp() {
    int width, height, channels;

    unsigned char *sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/openmp/Original.png", &width, &height, &channels,
                                     0);
    unsigned char *Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/openmp/Original.png", &width, &height,
                                            &channels, 0);
    unsigned char *Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/openmp/Original.png", &width, &height,
                                              &channels, 0);

    clock_t start = clock();

#pragma omp parallel for
    for (int i = 0; i < width * height * channels; i++) {
        Solarization[i] = solarize(Solarization[i]);
    }

#pragma omp parallel for
    for (int i = 0; i < width * height * channels; i++) {
        Pasteurization[i] = posterize(Pasteurization[i]);
    }

#pragma omp parallel for
    for (int i = 0; i < width * height * channels - 3; i = i + 3) {
        vector<unsigned char> rgb = toSepia(sepia[i], sepia[i + 1], sepia[i + 2]);
        sepia[i] = rgb[0];
        sepia[i + 1] = rgb[1];
        sepia[i + 2] = rgb[2];
    }

    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/openmp/Sepia.png", width, height, channels, sepia, 0)) {
        std::cout << "Sepia error" << std::endl;
        stbi_image_free(sepia);
    }
    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/openmp/Pasteurization.png", width, height, channels,
                        Pasteurization, 0)) {
        std::cout << "Pasteurization error" << std::endl;
        stbi_image_free(Pasteurization);
    }
    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/openmp/Solarization.png", width, height, channels,
                        Solarization, 0)) {
        
        stbi_image_free(Solarization);
    }
    return seconds;
}

double vectorize() {
    int width, height, channels;

    unsigned char *sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/vectorize/Original.png", &width, &height,
                                     &channels, 0);
    unsigned char *Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/vectorize/Original.png", &width, &height,
                                            &channels, 0);
    unsigned char *Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/vectorize/Original.png",
                                              &width, &height, &channels, 0);


    clock_t start = clock();

    //  Solarization
    __m128i threshold16 = _mm_set1_epi16(128);
    __m128i full255_16  = _mm_set1_epi16(255);

    for (int i = 0; i < width * channels * height / 8 * 8 - 8; i += 8) {
        __m128i bytes = _mm_loadl_epi64((__m128i const*)&Solarization[i]);
        __m128i vec16 = _mm_cvtepu8_epi16(bytes);

        __m128i mask = _mm_cmpgt_epi16(vec16, threshold16);

        __m128i subtracted = _mm_sub_epi16(full255_16, vec16);

        __m128i blended16 = _mm_blendv_epi8(vec16, subtracted, mask);

        __m128i final8 = _mm_packus_epi16(blended16, _mm_setzero_si128());

        _mm_storel_epi64((__m128i*)&Solarization[i], final8);
    }

    float step = 255.0f / 4.0f;

    int i = 0;
    for (; i + 8 < width * channels * height; i += 8) {

        __m128i pixels = _mm_loadl_epi64((__m128i*)(Pasteurization + i));

        __m128i pixels_16bit = _mm_cvtepu8_epi16(pixels);
        __m128i pixels_32bit_low = _mm_cvtepi16_epi32(pixels_16bit);
        __m128i pixels_32bit_high = _mm_cvtepi16_epi32(_mm_srli_si128(pixels_16bit, 8));

        __m128 pixels_float_low = _mm_cvtepi32_ps(pixels_32bit_low);
        __m128 pixels_float_high = _mm_cvtepi32_ps(pixels_32bit_high);

        __m128 step_vec = _mm_set1_ps(step);


        __m128 div_low = _mm_div_ps(pixels_float_low, step_vec);
        __m128 floor_low = _mm_floor_ps(div_low);
        __m128 result_low = _mm_mul_ps(floor_low, step_vec);

        __m128 div_high = _mm_div_ps(pixels_float_high, step_vec);
        __m128 floor_high = _mm_floor_ps(div_high);
        __m128 result_high = _mm_mul_ps(floor_high, step_vec);

        __m128i result_int_low = _mm_cvtps_epi32(result_low);
        __m128i result_int_high = _mm_cvtps_epi32(result_high);

        __m128i result_16bit = _mm_packus_epi32(result_int_low, result_int_high);

        __m128i result_8bit = _mm_packus_epi16(result_16bit, _mm_setzero_si128());

        _mm_storel_epi64((__m128i*)(Pasteurization + i), result_8bit);
    }

    for (; i < channels * width * height; i++) {
        Pasteurization[i] = (unsigned char)(int(int(Pasteurization[i]) / step) * step);
    }

    //

    const int pixels = width * height;
    i = 0;

    const __m128 r_coef = _mm_setr_ps(0.393f, 0.349f, 0.272f, 0.0f);
    const __m128 g_coef = _mm_setr_ps(0.769f, 0.686f, 0.534f, 0.0f);
    const __m128 b_coef = _mm_setr_ps(0.189f, 0.168f, 0.131f, 0.0f);
    const __m128 max_val = _mm_set1_ps(255.0f);


    for (; i < pixels; i++) {
        int idx = i * 3;

        unsigned char R = sepia[idx];
        unsigned char G = sepia[idx + 1];
        unsigned char B = sepia[idx + 2];

        __m128 r_vec = _mm_set1_ps((float)R);
        __m128 g_vec = _mm_set1_ps((float)G);
        __m128 b_vec = _mm_set1_ps((float)B);

        __m128 r_result = _mm_mul_ps(r_vec, r_coef);
        __m128 g_result = _mm_mul_ps(g_vec, g_coef);
        __m128 b_result = _mm_mul_ps(b_vec, b_coef);

        __m128 sum = _mm_add_ps(_mm_add_ps(r_result, g_result), b_result);

        __m128 clamped = _mm_min_ps(sum, max_val);

        float results[4];
        _mm_storeu_ps(results, clamped);

        sepia[idx]     = (unsigned char)results[0];
        sepia[idx + 1] = (unsigned char)results[1];
        sepia[idx + 2] = (unsigned char)results[2];
    }
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;

    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/vectorize/Solarization.png", width, height, channels,
                        Solarization, 0)) {
        stbi_image_free(Solarization);
    }

    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/vectorize/Pasteurization.png", width, height, channels,
                        Pasteurization, 0)) {
        stbi_image_free(Pasteurization);
                        }

    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/vectorize/Sepia.png", width, height, channels,
                        sepia, 0)) {
        stbi_image_free(sepia);
                        }

    return seconds;

}

int main() {

    int repeat = 1;

    double time1 = 0;
    for(int i = 0; i < repeat; i++) {
        time1 += openMp();
    }
    printf("openMp: %f seconds\n", time1);


    double time2 = 0;
    for(int i = 0; i < repeat; i++) {
        time2 += sequentially();
    }
    printf("sequentially: %f seconds\n", time2);

    double time3 = 0;
    for(int i = 0; i < repeat; i++) {
        time3 += vectorize();
    }
    printf("vectorize: %f seconds\n", time3);


}
