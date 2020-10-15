#include "helpers.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Convert image to sepia scale
void sepia(int height, int width, RGBTRIPLE image[height][width], RGBTRIPLE tmp_image[height][width], int y, int x)
{
  float sepiaRed = 0.393 * image[y][x].rgbtRed + 0.769 * image[y][x].rgbtGreen + 0.189 * image[y][x].rgbtBlue;
  float sepiaGreen = 0.349 * image[y][x].rgbtRed + 0.686 * image[y][x].rgbtGreen + 0.168 * image[y][x].rgbtBlue;
  float sepiaBlue = 0.272 * image[y][x].rgbtRed + 0.534 * image[y][x].rgbtGreen + 0.131 * image[y][x].rgbtBlue;

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

  tmp_image[y][x].rgbtRed = round(sepiaRed);
  tmp_image[y][x].rgbtGreen = round(sepiaGreen);
  tmp_image[y][x].rgbtBlue = round(sepiaBlue);
}

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width], RGBTRIPLE tmp_image[height][width], int y, int x)
{
  // getting average between RGB channels
  float avrg = (image[y][x].rgbtRed + image[y][x].rgbtGreen + image[y][x].rgbtBlue) / 3;

  // Assigning average value to all RGB channels
  tmp_image[y][x].rgbtRed = round(avrg);
  tmp_image[y][x].rgbtGreen = round(avrg);
  tmp_image[y][x].rgbtBlue = round(avrg);
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
  RGBTRIPLE tmp_image[height][width];

  // generating image copy
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      tmp_image[y][x] = image[y][x];
    }
  }

  int axis = 2;
  printf("============REFLECT OPTIONS=============\n");
  while (axis != 0 && axis != 1)
  {
    printf("Axis of reflection: (0 = X / 1 = Y): ");
    scanf("%d", &axis);
  }

  if (axis == 0)
  {
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

  printf("\nLoading...\n");

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
    // credit: https://www.geeksforgeeks.org/gaussian-filter-generation-c/
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
  return 0;
}

// Blurs an image
void blur(int height, int width, RGBTRIPLE image[height][width], RGBTRIPLE tmp_image[height][width], int y, int x, int kernel_size, double kernel[kernel_size][kernel_size])
{
  int kernel_limit = (kernel_size - 1) / 2;
  int opposite_kernel_limit = ((kernel_size - 1) / 2) * -1;

  float sumRed = 0;
  float sumGreen = 0;
  float sumBlue = 0;
  int neighbors = 0;

  for (int y_kernel = 0, y_kernel_aux = opposite_kernel_limit; y_kernel_aux <= kernel_limit; y_kernel++, y_kernel_aux++)
  {
    for (int x_kernel = 0, x_kernel_aux = opposite_kernel_limit; x_kernel_aux <= kernel_limit; x_kernel++, x_kernel_aux++)
    {
      int x_offset = x + x_kernel_aux;
      int y_offset = y + y_kernel_aux;

      if (x_offset < 0 || y_offset < 0 || x_offset >= width || y_offset >= height)
      {
        continue;
      }
      else
      {
        sumRed += image[y_offset][x_offset].rgbtRed * kernel[y_kernel][x_kernel];
        sumGreen += image[y_offset][x_offset].rgbtGreen * kernel[y_kernel][x_kernel];
        sumBlue += image[y_offset][x_offset].rgbtBlue * kernel[y_kernel][x_kernel];
        neighbors++;
      }
    }
  }
  int rgbtRed = round(sumRed / neighbors);
  int rgbtGreen = round(sumGreen / neighbors);
  int rgbtBlue = round(sumBlue / neighbors);

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

// Exposes edges of an image
void edge(int height, int width, RGBTRIPLE image[height][width], RGBTRIPLE tmp_image[height][width], int y, int x)
{
  int kernel_size = 3;
  // Sobel Operator X (Gx)
  int x_kernel_real[3][3] = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};

  // Sobel Operator Y (Gy)
  int y_kernel_real[3][3] = {{-1, -2, -1},
                             {0, 0, 0},
                             {1, 2, 1}};

  int kernel_limit = (kernel_size - 1) / 2;
  int opposite_kernel_limit = ((kernel_size - 1) / 2) * -1;

  float x_sumRed = 0;
  float x_sumGreen = 0;
  float x_sumBlue = 0;

  float y_sumRed = 0;
  float y_sumGreen = 0;
  float y_sumBlue = 0;

  int neighbors = 0;

  for (int y_kernel = 0, y_kernel_aux = opposite_kernel_limit; y_kernel_aux <= kernel_limit; y_kernel++, y_kernel_aux++)
  {
    for (int x_kernel = 0, x_kernel_aux = opposite_kernel_limit; x_kernel_aux <= kernel_limit; x_kernel++, x_kernel_aux++)
    {
      int x_offset = x + x_kernel_aux;
      int y_offset = y + y_kernel_aux;

      if (x_offset < 0 || y_offset < 0 || x_offset >= width || y_offset >= height)
      {
        continue;
      }
      else
      {
        // Gx
        x_sumRed += image[y_offset][x_offset].rgbtRed * x_kernel_real[y_kernel][x_kernel];
        x_sumGreen += image[y_offset][x_offset].rgbtGreen * x_kernel_real[y_kernel][x_kernel];
        x_sumBlue += image[y_offset][x_offset].rgbtBlue * x_kernel_real[y_kernel][x_kernel];

        // Gy
        y_sumRed += image[y_offset][x_offset].rgbtRed * y_kernel_real[y_kernel][x_kernel];
        y_sumGreen += image[y_offset][x_offset].rgbtGreen * y_kernel_real[y_kernel][x_kernel];
        y_sumBlue += image[y_offset][x_offset].rgbtBlue * y_kernel_real[y_kernel][x_kernel];
        neighbors++;
      }
    }
  }
  // Averaging neighbor pixels
  int gx_rgbtRed = round(x_sumRed / (float)neighbors);
  int gx_rgbtGreen = round(x_sumGreen / (float)neighbors);
  int gx_rgbtBlue = round(x_sumBlue / (float)neighbors);

  int gy_rgbtRed = round(y_sumRed / (float)neighbors);
  int gy_rgbtGreen = round(y_sumGreen / (float)neighbors);
  int gy_rgbtBlue = round(y_sumBlue / (float)neighbors);

  if (gx_rgbtRed > 255)
  {
    gx_rgbtRed = 255;
  }
  if (gx_rgbtGreen > 255)
  {
    gx_rgbtGreen = 255;
  }
  if (gx_rgbtBlue > 255)
  {
    gx_rgbtBlue = 255;
  }

  if (gy_rgbtRed > 255)
  {
    gy_rgbtRed = 255;
  }
  if (gy_rgbtGreen > 255)
  {
    gy_rgbtGreen = 255;
  }
  if (gy_rgbtBlue > 255)
  {
    gy_rgbtBlue = 255;
  }

  if (gx_rgbtRed < 0)
  {
    gx_rgbtRed *= -1;
  }
  if (gx_rgbtGreen < 0)
  {
    gx_rgbtGreen *= -1;
  }
  if (gx_rgbtBlue < 0)
  {
    gx_rgbtBlue *= -1;
  }

  if (gy_rgbtRed < 0)
  {
    gy_rgbtRed *= -1;
  }
  if (gy_rgbtGreen < 0)
  {
    gy_rgbtGreen *= -1;
  }
  if (gy_rgbtBlue < 0)
  {
    gy_rgbtBlue *= -1;
  }

  // G = (Gx^2 + Gy^2)^1/2
  tmp_image[y][x].rgbtRed = round(sqrt((gx_rgbtRed * gx_rgbtRed) + (gy_rgbtRed * gy_rgbtRed)));
  tmp_image[y][x].rgbtGreen = round(sqrt((gx_rgbtGreen * gx_rgbtGreen) + (gy_rgbtGreen * gy_rgbtGreen)));
  tmp_image[y][x].rgbtBlue = round(sqrt((gx_rgbtBlue * gx_rgbtBlue) + (gy_rgbtBlue * gy_rgbtBlue)));
}

// Final challenge!
void challenge(int height, int width, RGBTRIPLE image[height][width])
{
  RGBTRIPLE tmp_image[height][width];

  int y_quarter_limit = round(height / 2);

  reflect(height, width, image);

  // Getting the kernel_size correctly
  int kernel_size = 2;
  printf("\n============BLUR OPTIONS=============\n");
  while (kernel_size % 2 == 0 || kernel_size == 2)
  {
    printf("Kernel size: (only odd numbers): ");
    scanf("%d", &kernel_size);
  }

  double kernel[kernel_size][kernel_size];
  generate_blur_kernel(kernel_size, kernel);

  // Loop eixo y pela metade
  for (int y = 0; y < y_quarter_limit; y++)
  {
    // Blur completo
    for (int x = 0; x < width; x++)
    {
      blur(height, width, image, tmp_image, y, x, kernel_size, kernel);
    }
  }
// Loop eixo y pela outra metade
  for (int y = y_quarter_limit; y < height; y++)
  {
    // Sobel completo
    for (int x = 0; x < width; x++)
    {
      edge(height, width, image, tmp_image, y, x);
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