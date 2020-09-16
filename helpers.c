#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width])
{
    // height loop
    for (int i = 0; i < height; i++)
    {
        // width loop
        for (int j = 0; j < width; j++)
        {
            float avrg = (image[i][j].rgbtRed + image[i][j].rgbtGreen + image[i][j].rgbtBlue) / 3.0;

            int average = round(avrg);

            // Assigning the average between all RGB on the pixel
            image[i][j].rgbtBlue = average;
            image[i][j].rgbtRed = average;
            image[i][j].rgbtGreen = average;
        }
    }
}

// Convert image to sepia
void sepia(int height, int width, RGBTRIPLE image[height][width])
{
    // height loop
    for (int i = 0; i < height; i++)
    {
        // width loop
        for (int j = 0; j < width; j++)
        {
            float sepiaRed = (0.393 * image[i][j].rgbtRed) + (0.769 * image[i][j].rgbtGreen) + (0.189 * image[i][j].rgbtBlue);
            float sepiaGreen = (0.349 * image[i][j].rgbtRed) + (0.686 * image[i][j].rgbtGreen) + (0.168 * image[i][j].rgbtBlue);
            float sepiaBlue = (0.272 * image[i][j].rgbtRed) + (0.534 * image[i][j].rgbtGreen) + (0.131 * image[i][j].rgbtBlue);

            if (sepiaRed > 255)
            {
                sepiaRed = 255;
            }
            if (sepiaGreen > 255)
            {
                sepiaGreen = 255;
            }
            if (sepiaBlue > 255)
            {
                sepiaBlue = 255;
            }

            image[i][j].rgbtRed = round(sepiaRed);
            image[i][j].rgbtGreen = round(sepiaGreen);
            image[i][j].rgbtBlue = round(sepiaBlue);
        }
    }
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
    // used for 1. method
    RGBTRIPLE tmp_rgbt;

    // used in 2. method
    RGBTRIPLE tmp_image[height][width];

    // generating image copy
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            tmp_image[y][x] = image[y][x];
        }
    }

    int axis = 2;
    while (axis != 0 && axis != 1)
    {
        printf("Axis of reflection: (0 = X / 1 = Y): ");
        scanf("%d", &axis);
    }

    if (axis == 0)
    {
        // 1. METHOD
        // height loop
        // for (int y = 0; y < height; y++)
        // {
        //     // half width loop
        //     for (int x = 0; x < floor(width / 2); x++)
        //     {
        //         int reflected_index = (width - 1) - x;

        //         tmp_rgbt = image[y][reflected_index];
        //         image[y][reflected_index] = image[y][x];
        //         image[y][x] = tmp_rgbt;
        //     }
        // }

        // 2. METHOD
        // height loop
        for (int y = 0; y < height; y++)
        {
            // half width loop
            for (int x = 0; x < width; x++)
            {
                int reflected_index = (width - 1) - x;

                image[y][x] = tmp_image[y][reflected_index];
            }
        }
    }
    else if (axis == 1)
    {
        // 2. METHOD
        // height loop
        for (int y = 0; y < height; y++)
        {
            // half width loop
            for (int x = 0; x < width; x++)
            {
                int reflected_index = (height - 1) - y;

                image[y][x] = tmp_image[reflected_index][x];
            }
        }
    }
}

// generates kernel for mean and gaussian blur(0 = mean blur / 1 = gaussian blur)
int generate_blur_kernel(size_t kernel_size, double kernel[kernel_size][kernel_size])
{
    // Getting type correctly
    int type = 2;
    while (type != 0 && type != 1)
    {
        printf("Type of kernel: (0 = Mean blur / 1 = Gaussian blur): ");
        scanf("%d", &type);
    }

    if (type == 0)
    {
        // Mean blur (filling all the positions with 1)
        for (int y = 0; y < kernel_size; y++)
        {
            for (int x = 0; x < kernel_size; x++)
            {
                kernel[y][x] = 1;
            }
        }
    }
    else if (type == 1)
    {
        int kernel_iterator = (kernel_size - 1) / 2;
        int opposite_kernel_iterator = ((kernel_size - 1) / 2) * -1;

        // Gaussian filter kernel generator
        // (most haha) credit: https://www.geeksforgeeks.org/gaussian-filter-generation-c/
        double sigma = 1.0;
        double r, s = 2.0 * sigma * sigma;
        #define M_PI 3.14159265358979323846

        // sum for normalization
        double sum = 0.0;

        for (int x = opposite_kernel_iterator; x <= kernel_iterator; x++)
        {
            for (int y = opposite_kernel_iterator; y <= kernel_iterator; y++)
            {
                r = sqrt(x * x + y * y);
                kernel[x + kernel_iterator][y + kernel_iterator] = (exp(-(r * r) / s)) / (M_PI * s);
                sum += kernel[x + kernel_iterator][y + kernel_iterator];
            }
        }

        // normalizing the Kernel
        for (int i = 0; i < kernel_size; ++i)
        {
            for (int j = 0; j < kernel_size; ++j)
            {
                kernel[i][j] = (kernel[i][j] / sum) + 1;
            }
        }
    }

    // Printing the kernel
    printf("Your kernel:\n");
    for (int y = 0; y < kernel_size; y++)
    {
        for (int x = 0; x < kernel_size; x++)
        {
            printf("%f  ", kernel[y][x]);
        }
        printf("\n");
    }
    return 0;
}

// Blur image
void blur(int height, int width, RGBTRIPLE image[height][width])
{
    RGBTRIPLE tmp_image[height][width];

    // Getting the kernel_size correctly
    static int kernel_size;
    while (kernel_size % 2 == 0 && kernel_size == 0)
    {
        printf("Kernel size: (only odd numbers): ");
        scanf("%d", &kernel_size);
    }

    double kernel[kernel_size][kernel_size];
    generate_blur_kernel(kernel_size, kernel);

    int kernel_iterator = (kernel_size - 1) / 2;
    int opposite_kernel_iterator = ((kernel_size - 1) / 2) * -1;

    // heigth loop
    for (int y = 0; y < height; y++)
    {
        // width loop
        for (int x = 0; x < width; x++)
        {
            // Our own typedef for operations
            float_RGBTRIPLE tmp_avrg;
            // Initializing our typedef preveting from errors
            tmp_avrg.rgbtRed = 0;
            tmp_avrg.rgbtGreen = 0;
            tmp_avrg.rgbtBlue = 0;

            // Neighbors count to dividing in the end
            int neighbor_count = 0;

            // Neighbors sum
            for (int o = opposite_kernel_iterator, o_kernel = 0; o <= kernel_iterator; o++, o_kernel++)
            {
                for (int u = opposite_kernel_iterator, u_kernel = 0; u <= kernel_iterator; u++, u_kernel++)
                {
                    int x_offset = x + u;
                    int y_offset = y + o;

                    if (x_offset < 0 || x_offset >= width || y_offset < 0 || y_offset >= height)
                    {
                        continue;
                    }
                    else
                    {
                        tmp_avrg.rgbtRed += image[y_offset][x_offset].rgbtRed * kernel[o_kernel][u_kernel];
                        tmp_avrg.rgbtGreen += image[y_offset][x_offset].rgbtGreen * kernel[o_kernel][u_kernel];
                        tmp_avrg.rgbtBlue += image[y_offset][x_offset].rgbtBlue * kernel[o_kernel][u_kernel];
                        neighbor_count++;
                    }
                }
            }

            int rgbtRed = round(tmp_avrg.rgbtRed / (float)neighbor_count);
            int rgbtGreen = round(tmp_avrg.rgbtGreen / (float)neighbor_count);
            int rgbtBlue = round(tmp_avrg.rgbtBlue / (float)neighbor_count);

            if (rgbtRed > 255)
            {
                rgbtRed = 255;
            }
            if (rgbtGreen > 255)
            {
                rgbtGreen = 255;
            }
            if (rgbtBlue > 255)
            {
                rgbtBlue = 255;
            }

            tmp_image[y][x].rgbtRed = rgbtRed;
            tmp_image[y][x].rgbtGreen = rgbtGreen;
            tmp_image[y][x].rgbtBlue = rgbtBlue;
        }
    }

    // copying pixels into original image
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            image[y][x] = tmp_image[y][x];
        }
    }
}

// Exposes edges of an image
void edge(int height, int width, RGBTRIPLE image[height][width])
{
    RGBTRIPLE tmp_image[height][width];

    // Sobel Operator X (Gx)
    int x_kernel[3][3] = {{-1, 0, 1},
                          {-2, 0, 2},
                          {-1, 0, 1}};

    // Soble Operator Y (Gy)
    int y_kernel[3][3] = {{-1, -2, -1},
                          {0, 0, 0},
                          {1, 2, 1}};

    // heigth loop
    for (int y = 0; y < height; y++)
    {
        // width loop
        for (int x = 0; x < width; x++)
        {

            float x_tmp_avrgRed = 0;
            float x_tmp_avrgGreen = 0;
            float x_tmp_avrgBlue = 0;

            float y_tmp_avrgRed = 0;
            float y_tmp_avrgGreen = 0;
            float y_tmp_avrgBlue = 0;

            int neighbor_count = 0;

            // Neighbors sum
            for (int o = -1, o_kernel = 0; o <= 1; o++, o_kernel++)
            {
                for (int u = -1, u_kernel = 0; u <= 1; u++, u_kernel++)
                {
                    int x_offset = x + u;
                    int y_offset = y + o;

                    if (x_offset < 0 || x_offset >= width || y_offset < 0 || y_offset >= height)
                    {
                        continue;
                    }
                    else
                    {
                        // Gy
                        x_tmp_avrgRed += image[y_offset][x_offset].rgbtRed * x_kernel[o_kernel][u_kernel];
                        x_tmp_avrgGreen += image[y_offset][x_offset].rgbtGreen * x_kernel[o_kernel][u_kernel];
                        x_tmp_avrgBlue += image[y_offset][x_offset].rgbtBlue * x_kernel[o_kernel][u_kernel];

                        // Gx
                        y_tmp_avrgRed += image[y_offset][x_offset].rgbtRed * y_kernel[o_kernel][u_kernel];
                        y_tmp_avrgGreen += image[y_offset][x_offset].rgbtGreen * y_kernel[o_kernel][u_kernel];
                        y_tmp_avrgBlue += image[y_offset][x_offset].rgbtBlue * y_kernel[o_kernel][u_kernel];

                        neighbor_count++;
                    }
                }
            }

            // Averaging neighbor pixels
            int x_rgbtRed = round(x_tmp_avrgRed / (float)neighbor_count);
            int x_rgbtGreen = round(x_tmp_avrgGreen / (float)neighbor_count);
            int x_rgbtBlue = round(x_tmp_avrgBlue / (float)neighbor_count);

            int y_rgbtRed = round(y_tmp_avrgRed / (float)neighbor_count);
            int y_rgbtGreen = round(y_tmp_avrgGreen / (float)neighbor_count);
            int y_rgbtBlue = round(y_tmp_avrgBlue / (float)neighbor_count);

            if (x_rgbtRed > 255)
            {
                x_rgbtRed = 255;
            }
            if (x_rgbtGreen > 255)
            {
                x_rgbtGreen = 255;
            }
            if (x_rgbtBlue > 255)
            {
                x_rgbtBlue = 255;
            }

            if (y_rgbtRed > 255)
            {
                y_rgbtRed = 255;
            }
            if (y_rgbtGreen > 255)
            {
                y_rgbtGreen = 255;
            }
            if (y_rgbtBlue > 255)
            {
                y_rgbtBlue = 255;
            }

            if (x_rgbtRed < 0)
            {
                x_rgbtRed *= -1;
            }
            if (x_rgbtGreen < 0)
            {
                x_rgbtGreen *= -1;
            }
            if (x_rgbtBlue < 0)
            {
                x_rgbtBlue *= -1;
            }

            if (y_rgbtRed < 0)
            {
                y_rgbtRed *= -1;
            }
            if (y_rgbtGreen < 0)
            {
                y_rgbtGreen *= -1;
            }
            if (y_rgbtBlue < 0)
            {
                y_rgbtBlue *= -1;
            }

            // G = (Gx^2 + Gy^2)^1/2
            // Convolution
            int rgbtRed = sqrt((x_rgbtRed * x_rgbtRed) + (y_rgbtRed * y_rgbtRed));
            int rgbtGreen = sqrt((x_rgbtGreen * x_rgbtGreen) + (y_rgbtGreen * y_rgbtGreen));
            int rgbtBlue = sqrt((x_rgbtBlue * x_rgbtBlue) + (y_rgbtBlue * y_rgbtBlue));

            tmp_image[y][x].rgbtRed = round(rgbtRed);
            tmp_image[y][x].rgbtGreen = round(rgbtGreen);
            tmp_image[y][x].rgbtBlue = round(rgbtBlue);
        }
    }

    // copying pixels into original image
    /*
    1. We can optmize this part to "cut" the white borders of the image,
    but we are not doing this here because of the architecture of filter.c

    2. Another solution would be create a "padding" (extend) for our picture,
    this way we recreate a small (1, 2, 3 px) border with the most possible 
    colors for the border.
    We're using this tchnique here ^^^
    */
    for (int y = 0; y < height; y++)
    {
        int y_offset = y;
        if (y == 0)
        {
            y_offset += 1;
        }
        else if (y == height - 1)
        {
            y_offset -= 1;
        }

        for (int x = 0; x < width; x++)
        {
            int x_offset = x;
            if (x == 0)
            {
                x_offset += 1;
            }
            else if (x == width - 1)
            {
                x_offset -= 1;
            }

            image[y][x] = tmp_image[y_offset][x_offset];
        }
    }
}