#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <time.h>

#define SHIFT 2 // if 5x5 kernel; 1, if 3x3 kernel

void convert_to_gray(unsigned char* image, unsigned width, unsigned height) {
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x) {
            unsigned char r = image[4 * width * y + 4 * x + 0];
            unsigned char g = image[4 * width * y + 4 * x + 1];
            unsigned char b = image[4 * width * y + 4 * x + 2];
            unsigned char gray = (unsigned char)(0.21 * r + 0.72 * g + 0.07 * b);
            image[4 * width * y + 4 * x + 0] = gray;
            image[4 * width * y + 4 * x + 1] = gray;
            image[4 * width * y + 4 * x + 2] = gray;
            image[4 * width * y + 4 * x + 2] = 255;
        }
    }
}

void gaussian_filter(unsigned char* image, unsigned width, unsigned height) {
    unsigned char* temp_image = malloc(4 * width * height * sizeof(unsigned char));
    if (!temp_image) {
        fprintf(stderr, "Error: Memory allocation failed for temp_image\n");
        exit(1);
    }
    float kernel[3][3] = {
        {1.0 / 16, 2.0 / 16, 1.0 / 16},
        {2.0 / 16, 4.0 / 16, 2.0 / 16},
        {1.0 / 16, 2.0 / 16, 1.0 / 16}
    };
    for (unsigned y = 1; y < height - 1; ++y) {
        for (unsigned x = 1; x < width - 1; ++x) {
            float sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    unsigned char r = image[4 * width * (y + i) + 4 * (x + j) + 0];
                    unsigned char g = image[4 * width * (y + i) + 4 * (x + j) + 1];
                    unsigned char b = image[4 * width * (y + i) + 4 * (x + j) + 2];
                    sum_r += kernel[i + 1][j + 1] * r;
                    sum_g += kernel[i + 1][j + 1] * g;
                    sum_b += kernel[i + 1][j + 1] * b;
                }
            }
            temp_image[4 * width * y + 4 * x + 0] = (unsigned char)fmin(sum_r, 255);
            temp_image[4 * width * y + 4 * x + 1] = (unsigned char)fmin(sum_g, 255);
            temp_image[4 * width * y + 4 * x + 2] = (unsigned char)fmin(sum_b, 255);
        }
    }
    for (unsigned y = 1; y < height - 1; ++y) {
        for (unsigned x = 1; x < width - 1; ++x) {
            image[4 * width * y + 4 * x + 0] = temp_image[4 * width * y + 4 * x + 0];
            image[4 * width * y + 4 * x + 1] = temp_image[4 * width * y + 4 * x + 1];
            image[4 * width * y + 4 * x + 2] = temp_image[4 * width * y + 4 * x + 2];
        }
    }
    free(temp_image);
}


void sobel_operator(unsigned char* image, unsigned width, unsigned height) {
    unsigned char* temp_image = malloc(4 * width * height * sizeof(unsigned char));
    if (!temp_image) {
        fprintf(stderr, "Error: Memory allocation failed for temp_image\n");
        exit(1);
    }
    /*int sobel_x[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    int sobel_y[3][3] = {
        {-1, -2, -1},
        { 0,  0,  0},
        { 1,  2,  1}
    };*/
    const char sobel_x[5][5] = {
        {+2, +1, 0, -1, -2},
        {+2, +1, 0, -1, -2},
        {+4, +2, 0, -2, -4},
        {+2, +1, 0, -1, -2},
        {+2, +1, 0, -1, -2}
    };
    const char sobel_y[5][5] = {
        {+2, +2, +4, +2, +2},
        {+1, +1, +2, +1, +1},
        {+0, +0, +0, +0, +0},
        {-1, -1, -2, -1, -1},
        {-2, -2, -4, -2, -2}
    };
    for (unsigned y = SHIFT; y < height - SHIFT; ++y) {
        for (unsigned x = SHIFT; x < width - SHIFT; ++x) {
            float sum_r_x = 0.0, sum_g_x = 0.0, sum_b_x = 0.0;
            float sum_r_y = 0.0, sum_g_y = 0.0, sum_b_y = 0.0;
            for (int i = -SHIFT; i <= SHIFT; ++i) {
                for (int j = -SHIFT; j <= SHIFT; ++j) {
                    unsigned char r = image[4 * width * (y + i) + 4 * (x + j) + 0];
                    unsigned char g = image[4 * width * (y + i) + 4 * (x + j) + 1];
                    unsigned char b = image[4 * width * (y + i) + 4 * (x + j) + 2];
                    sum_r_x += sobel_x[i + SHIFT][j + SHIFT] * r;
                    sum_g_x += sobel_x[i + SHIFT][j + SHIFT] * g;
                    sum_b_x += sobel_x[i + SHIFT][j + SHIFT] * b;
                    sum_r_y += sobel_y[i + SHIFT][j + SHIFT] * r;
                    sum_g_y += sobel_y[i + SHIFT][j + SHIFT] * g;
                    sum_b_y += sobel_y[i + SHIFT][j + SHIFT] * b;
                }
            }
            float magnitude_r = sqrt(sum_r_x * sum_r_x + sum_r_y * sum_r_y);
            float magnitude_g = sqrt(sum_g_x * sum_g_x + sum_g_y * sum_g_y);
            float magnitude_b = sqrt(sum_b_x * sum_b_x + sum_b_y * sum_b_y);
            temp_image[4 * width * y + 4 * x + 0] = (unsigned char)fmin(magnitude_r, 255);
            temp_image[4 * width * y + 4 * x + 1] = (unsigned char)fmin(magnitude_g, 255);
            temp_image[4 * width * y + 4 * x + 2] = (unsigned char)fmin(magnitude_b, 255);
        }
    }
    for (unsigned y = SHIFT; y < height - SHIFT; ++y) {
        for (unsigned x = SHIFT; x < width - SHIFT; ++x) {
            image[4 * width * y + 4 * x + 0] = temp_image[4 * width * y + 4 * x + 0];
            image[4 * width * y + 4 * x + 1] = temp_image[4 * width * y + 4 * x + 1];
            image[4 * width * y + 4 * x + 2] = temp_image[4 * width * y + 4 * x + 2];
        }
    }
    free(temp_image);
}

typedef struct Coords {
    unsigned x, y;
} Coords;

bool valid_coords(unsigned x, unsigned y, unsigned width, unsigned height) {
    return x >= SHIFT && y >= SHIFT && x < width - SHIFT && y < height - SHIFT;
}

bool is_border(unsigned char* image, unsigned x, unsigned y, unsigned width) {
    return image[4 * (y * width + x)] > 150; // 20, if 3x3 kernel; >= 150, if 5x5 kernel
}

void bfs(unsigned char* output, unsigned char* image, unsigned width, unsigned height, unsigned x, unsigned y, bool* visited) {
    unsigned char red = 120 + abs(rand()) % 136;
    unsigned char green = 180 + abs(rand()) % 76;
    unsigned char blue = 120 + abs(rand()) % 136;
    Coords* queue = malloc(width * height * sizeof(Coords));
    unsigned left = 0, right = 1;
    Coords init = { x, y };
    queue[0] = init;
    while (right - left > 0) {
        Coords current = queue[left++];
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                unsigned newx = current.x + i, newy = current.y + j;
                if (!valid_coords(newx, newy, width, height)) {
                    continue;
                }
                if (visited[newy * width + newx]) {
                    continue;
                }
                if (is_border(image, newx, newy, width)) {
                    continue;
                }
                visited[newy * width + newx] = true;
                Coords newcoords = { newx, newy };
                queue[right++] = newcoords;
            }
        }
        output[4 * (current.y * width + current.x) + 0] = red;
        output[4 * (current.y * width + current.x) + 1] = green;
        output[4 * (current.y * width + current.x) + 2] = blue;
        output[4 * (current.y * width + current.x) + 3] = 255;
    }
    free(queue);
    return;
}

void segmentation(unsigned char* output, unsigned char* image, unsigned width, unsigned height) {
    bool* visited = calloc(width * height, sizeof(bool));
    for (unsigned y = SHIFT; y < height - SHIFT; ++y) {
        for (unsigned x = SHIFT; x < width - SHIFT; ++x) {
            if (is_border(image, x, y, width)) {
                output[4 * (y * width + x) + 0] = 0;
                output[4 * (y * width + x) + 1] = 0;
                output[4 * (y * width + x) + 2] = 0;
                output[4 * (y * width + x) + 3] = 255;
                continue;
            }
            if (visited[width * y + x]) {
                continue;
            }
            bfs(output, image, width, height, x, y, visited);
        }
    }
    free(visited);
    return;
}

int main() {
    srand(8);
    const char* input_filename = "head.png";
    const char* output_filename = "output_scull.png";
    unsigned char* segmented_image, * colored_image;
    unsigned error;
    unsigned char* image;
    unsigned width, height;

    error = lodepng_decode32_file(&image, &width, &height, input_filename);
    if (error) {
        printf("error %u: %s\n", error, lodepng_error_text(error));
        return 1;
    }
    unsigned char* output = calloc(4 * width * height, sizeof(unsigned char));
    convert_to_gray(image, width, height);
    gaussian_filter(image, width, height);
    sobel_operator(image, width, height);
    for (unsigned y = 1; y < height - 1; ++y) {
        for (unsigned x = 1; x < width - 1; ++x) {
            output[4 * (y * width + x) + 0] = image[4 * (y * width + x) + 0];
            output[4 * (y * width + x) + 1] = image[4 * (y * width + x) + 1];
            output[4 * (y * width + x) + 2] = image[4 * (y * width + x) + 2];
            output[4 * (y * width + x) + 3] = 255;
        }
    }
    error = lodepng_encode32_file("sobel.png", output, width, height);
    segmentation(output, image, width, height);
    error = lodepng_encode32_file(output_filename, output, width, height);
    if (error) {
        printf("error %u: %s\n", error, lodepng_error_text(error));
        return 1;
    }

    free(image);
    free(output);
    return 0;
}
