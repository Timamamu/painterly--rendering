

#include "Image.h"
#include "basicImageManipulation.h"
#include "npr.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;


void testPaint() {
  // load images
  Image input_image("./Input/villeperdue.png");
  Image texture_imaage("./Input/brush.png");

  //autostitchN(ims, 1).write("./Output/myimages-autostitchN.png");
  //brush(Image &im, int x, int y, vector<float> color, const Image &texture)

  //get dimensions 
  int width = input_image.width();
  int height = input_image.height();
  int channels = input_image.channels();

  vector<float> color = {0.0f, 0.0f, 1.0f};
  //Image output(width, height, channels);
  // Create an output image (copy of the input)
  Image output = input_image;

  //random number generator
  srand(static_cast<unsigned int>(time(0)));

  //apply 10 random brush strokes 
  for (int i = 0; i < 10; ++i) {
    int x = rand() % width;
    int y = rand() % height;

    brush(output, x, y, color, texture_imaage);
  }

  output.write("./Output/testPaint_output.png");

}


void testSingleScalePaint() {
    // Load input image and texture (brush shape)
    Image input_image("./Input/villeperdue.png");
    Image texture_image("./Input/brush.png");

    // Define parameters for the painting effect
    int size = 10;         
    int N = 10000;          
    float noise = 0.3f;    

    // Create an output image as a copy of the input image to apply the brush strokes onto
    Image output_image = input_image;

    // Apply the single-scale paint effect
    //singleScalePaint(const Image &im, Image &out, const Image &texture, int size, int N, float noise)
    //singleScalePaint()
    singleScalePaint(input_image, output_image, texture_image, size, N, noise);

    // Save the output image
    output_image.write("./Output/testSingleScalePaint_output.png");

    
}

void testSingleScalePaintImportance() {
    // Load the input image and brush texture
    Image input_image("./Input/villeperdue.png");   // Main image to paint over
    Image texture_image("./Input/brush.png");       // Brush texture shape

    // Generate or load an importance map
    // Here, we assume a grayscale image as an importance map with values between 0 and 1.
    Image importance_map(input_image.width(), input_image.height(), 1);
    
    // Fill the importance map with a gradient (for example, high importance on the left, low on the right)
    for (int y = 0; y < importance_map.height(); ++y) {
        for (int x = 0; x < importance_map.width(); ++x) {
            float importance_value = static_cast<float>(x) / importance_map.width(); // Left (high) to right (low)
            importance_map(x, y, 0) = importance_value;
        }
    }

    // Define parameters for the painting effect
    int size = 10;         
    int N = 10000;          
    float noise = 0.3f;   

    // Create an output image as a copy of the input image
    Image output_image = input_image;

    // Apply the single-scale paint effect with importance sampling
    singleScalePaintImportance(input_image, importance_map, output_image, texture_image, size, N, noise);

    // Save the output image
    output_image.write("./Output/testSingleScalePaintImportance_output.png");

    
}

void testpainterly() {
  Image input("./Input/villeperdue.png");  
  Image texture("./Input/brush.png");      
  Image output(input.width(), input.height(), input.channels());                    

  int N = 10000;          // Number of strokes
  int size = 50;         // Maximum size of the brush
  float noise = 0.3f;    // Noise factor for color variation

  painterly(input, output, texture, N, size, noise);
}

void testOrientationExtraction() {
    
    Image input("./Input/round.png");
    float sigmaG = 3.0f;    
    float factor = 5.0f;
    // float sigmaG = 1.0f;    
    // float factor = 4.0f;      
    Image orientationMap = testAngle(input, sigmaG, factor);
    orientationMap.write("./Output/orientation_map.png");
}

void testSingleScaleOrientedPaint() {
    Image input("./Input/china.png");
    Image texture("./Input/brush.png");

    // Parameters
    int size = 10;            
    int N = 10000;             
    float noise = 0.3f;       
    int nAngles = 16;         

    // Output image
    Image output(input.width(), input.height(), input.channels());

    //Compute importance map (constant for testing purposes)
    Image importance(input.width(), input.height(), 1);
    importance.set_color(1.0f);  

    // Compute structure tensor
    float sigmaG = 1.0f;
    float factorSigma = 4.0f;
    Image tensor = computeTensor(input, sigmaG, factorSigma);

    // Run single-scale oriented paint function
    singleScaleOrientedPaint(input, importance, output, tensor, texture, size, N, noise, nAngles);

    // Save output image
    output.write("./Output/oriented_paint_china.png");  
}

void testOrientedPaint() {
    
    Image input("./Input/my_image.png");
    Image texture("./Input/brush.png");
    Image output(input.width(), input.height(), input.channels());

    // Define parameters for the oriented painterly effect
    int N = 10000;        
    int size = 50;        
    float noise = 0.3f;   
    orientedPaint(input, output, texture, N, size, noise);
    //output.write("./Output/oriented_paint_test_output.png");

    
}

void testPaintTransformer() {
    
    Image target("./Input/paper_image.png");
    Image texture("./Input/brush.png");
    //Image output(input.width(), input.height(), input.channels());

    // Set up parameters for the paint transformer function
    int numStrokes = 1000;  // Adjust based on desired density
    float noiseLevel = 0.1f;
    int numScales = 3;

    Image canvas(target.width(), target.height(), target.channels());

    //paintTransformer(target, canvas, texture, numStrokes, noiseLevel, numScales);
    paintTransformer(target, canvas, numStrokes, noiseLevel, numScales);
    //canvas.write("./Output/paint_transformer_output.png");

    
}

int main() {
  

  //testPaint();
  //testSingleScalePaint();
  //testSingleScalePaintImportance();
  //testpainterly();
  //testOrientationExtraction();
  //testSingleScaleOrientedPaint();
  //testOrientedPaint();
  testPaintTransformer();

  return 0;
}
