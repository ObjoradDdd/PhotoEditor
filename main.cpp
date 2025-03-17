#define STB_IMAGE_IMPLEMENTATION
#include "libs/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "libs/stb_image_write.h"
#include <iostream>
#include <vector>

using namespace std;

unsigned char solarize(unsigned char value){
    if ((int)value > 128) {
        return  (unsigned char)(255 - (int)value);
    }
    return value;
}

unsigned char posterize(unsigned char value, int levels = 5) {
    float step = 255.0f / (levels - 1);
    return (unsigned char)(int(int(value) / step) * step);
}

vector<unsigned char> toSepia(unsigned char R,unsigned char G, unsigned char B) {
    vector<unsigned char> rgb(3);
    rgb[0] = (unsigned char)(int(R) * 0.393 + int(G) * 0.769 + int(B) * 0.189);
    rgb[1] = (unsigned char)(int(R) * 0.349 + int(G) * 0.686 + int(B) * 0.168);
    rgb[2] = (unsigned char)(int(R) * 0.272 + int(G) * 0.534 + int(B) * 0.131);

    return rgb;
}

void sequentially() {
    int width, height, channels;

    unsigned char* sepia = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width, &height, &channels, 0);
    unsigned char* Solarization = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width, &height, &channels, 0);
    unsigned char* Pasteurization = stbi_load("C:/Users/Atrun/CLionProjects/fe/sequentially/Original.png", &width, &height, &channels, 0);

    if (!sepia) {
        std::cout << "Downloading error" << std::endl;
        return;
    }

    for (int i = 0; i < width * height * channels; i++) {
        Solarization[i] = solarize(Solarization[i]);
    }

    for (int i = 0; i < width * height * channels; i++) {
        Pasteurization[i] = posterize(Pasteurization[i]);
    }

    for (int i = 0; i < width * height * channels - 3; i = i+3) {
        vector<unsigned char> rgb = toSepia(sepia[i], sepia[i+1], sepia[i+2]);
        sepia[i] = rgb[0];
        sepia[i+1] = rgb[1];
        sepia[i+2] = rgb[2];
    }


    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Sepia.png", width, height, channels, sepia, 0)) {
        std::cout << "Sepia error" << std::endl;
        stbi_image_free(sepia);
    }
    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Pasteurization.png", width, height, channels, Pasteurization, 0)) {
        std::cout << "Pasteurization error" << std::endl;
        stbi_image_free(Pasteurization);
    }
    if (!stbi_write_png("C:/Users/Atrun/CLionProjects/fe/sequentially/Solarization.png", width, height, channels, Solarization, 0)) {
        std::cout << "Solarization error" << std::endl;
        stbi_image_free(Solarization);
    }
    
    
}


int main() {
    sequentially();
}