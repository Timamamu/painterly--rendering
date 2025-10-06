// npr.cpp

#include "npr.h"
#include "filtering.h"
#include "matrix.h"
#include <algorithm>
#include <math.h>
#include <random>
//#include <cmath>

using namespace std;

/**************************************************************
 //                       NPR                                //
 *************************************************************/

void brush(Image &im, int x, int y, vector<float> color, const Image &texture) {
	// Draws a brushstroke defined by texture and color at (x,y) in im
	int halfWidth = texture.width() / 2;
	int halfHeight = texture.height() / 2;

	//check boundaries
	if(x < halfWidth || x >= im.width() - halfWidth || y < halfHeight || y >= im.height() - halfHeight) {
		return;
	}

	//apply brush stroke
	for (int j = 0; j < texture.height(); ++j) {
		for (int i = 0; i < texture.width(); ++i) {
				//compute position 
				int im_x = x - halfWidth + i;
				int im_y = y - halfHeight + j;

				//get opacity from texture
				float opacity = texture(i, j, 0);

				//blend each color channel of the image
				for (int c = 0; c < im.channels(); ++c) {
					im(im_x, im_y, c) = (1 - opacity) * im(im_x, im_y, c) + opacity * color[c];
				}
		}
	}
	
}

void singleScalePaint(const Image &im, Image &out, const Image &texture, int size, int N, float noise) {
	// Create painted rendering by splatting brushstrokes at N random locations
	// in your ouptut image
	float factor = (texture.width() > texture.height()) 
                   ? static_cast<float>(size) / texture.width() 
                   : static_cast<float>(size) / texture.height();


	//random number generator
	srand(static_cast<unsigned int>(time(0)));

	//scale the texture to have maximum size of size
	Image scaledTexture = scaleLin(texture, factor);

	int width = out.width();
	int height = out.height();

	for (int i = 0; i < N; ++i) {
		//generate random location
		int x = rand() % width;
		int y = rand() % height;

		//get the color at this random location from im 
		vector<float> color(im.channels());
		for (int c = 0; c < im.channels(); ++c) {
			color[c] = im(x, y, c);
		}

		//adjust color by a noise factor 
		float noiseFactor = 1 - noise / 2 + noise * ((float)rand() / RAND_MAX);
		for (int c = 0; c < color.size(); ++c) {
			color[c] *= noiseFactor;
		}

		//apply brush stroke at x,y 
		brush(out, x, y, color, scaledTexture);
	}
	
}

void singleScalePaintImportance(const Image &im, const Image &importance,
						Image &out, const Image &texture, int size, int N, float noise) {
	// Create painted rendering but vary the density of the strokes according to
	// an importance map
	float avgImportance = 0.0f;
	int width = importance.width();
	int height = importance.height();

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			avgImportance += importance(x, y, 0);
		}
	}
	avgImportance /= (width * height);

	//adjust N based on the average importance
	int adjustedN = static_cast<int> (N / avgImportance);

	//scale the texture to have a maximum dimension size
	float factor = (texture.width() > texture.height())
					? static_cast<float>(size) / texture.width()
					: static_cast<float>(size) / texture.height();

	Image scaledTexture = scaleLin(texture, factor);

	//random number generator
	srand(static_cast<unsigned int>(time(0)));

	int outWidth = out.width();
	int outHeight = out.height();

	int strokesApplied = 0;

	//apply brush strokes based on importance
	for (int i = 0; i < adjustedN; ++i) {
		//randomly choose a position
		int x = rand() % outWidth;
		int y = rand() % outHeight;

		//get the importance value at this location 
		float importanceValue = importance(x, y, 0);

		//accept sample with probability importance value
		float randomValue = static_cast<float>(rand()) / RAND_MAX;
		if (randomValue > importanceValue) {
			//reject sample
			continue;
		}

		//get color from input image
		vector<float>color(im.channels());
		for (int c = 0; c < im.channels(); ++c) {
			color[c] = im(x, y, c);
		}

		//adjust the color by a noise factor 
		float noiseFactor = 1 - noise / 2 + noise * ((float)rand() / RAND_MAX);
		for (int c = 0; c < color.size(); ++c) {
			color[c] *= noiseFactor;
		}

		//apply brush stroke
		brush(out, x, y, color, scaledTexture);

		//increment strokesapplied
		strokesApplied++;
		if (strokesApplied >= N) {
			break;
		}
	}

	return;
}

Image sharpnessMap(const Image &im, float sigma) {
	// Calculate sharpness mask 
	
	//convert image to luminance
	vector<Image>lumiChromV = lumiChromi(im);
	Image luminance = lumiChromV[0];

	//apply guassian blur to luminance
	Image blurredLumi = gaussianBlur_separable(luminance, sigma);

	//calculate high frequencies(luminance - blurred luminance)
	Image highFreq = luminance - blurredLumi;

	//sqaure the high freq image 
	for (int y = 0; y < highFreq.height(); ++y) {
		for (int x = 0; x < highFreq.width(); ++x) {
			highFreq(x, y, 0) = highFreq(x, y, 0) * highFreq(x, y, 0);
		}
	}

	//apply gussian blur to the squared high freq image with 4 * sigma
	Image blurredHighFreq = gaussianBlur_separable(highFreq, 4 * sigma);

	//normalise by the max value
	float maxVal = 0.0f;
	for (int y = 0; y < blurredHighFreq.height(); ++y) {
		for (int x = 0; x < blurredHighFreq.width(); ++x) {
			maxVal = max(maxVal, blurredHighFreq(x, y, 0));
		}
	}

	//normalize but ensure not by zero
	if (maxVal > 0) {
		for (int y = 0; y < blurredHighFreq.height(); ++y) {
			for (int x = 0; x < blurredHighFreq.width(); ++x) {
				blurredHighFreq(x, y, 0) /= maxVal;
			}
		}
	}
	return blurredHighFreq;
}

void painterly(const Image &im, Image &out, const Image &texture, int N, int size, float noise) {
	// Create painterly rendering using a first layer of coarse strokes followed
	// by smaller strokes in high detail areas
	
	//create a constant importance map for the first pass
	Image constantImportanceMap (im.width(), im.height(), 1);
	constantImportanceMap.set_color(1.0f);

	//apply first layer of coarse strokes with constant importance
	int coarseSize = size;
	singleScalePaintImportance(im, constantImportanceMap, out, texture, coarseSize, N, noise);

	//save after first pass
	out.write("./Output/painterly_output_first_pass.png");

	//generate sharpness map for the second pass
	float sigma = 1.0f;
	Image sharpness = sharpnessMap(im, sigma);
	

	//apply seond layer of finer strokes
	int fineSize = size / 4;
	int refinedN = N / 2;
	singleScalePaintImportance(im, sharpness, out, texture, fineSize, refinedN, noise);

	out.write("./Output/painterly_output_final.png");
	//return;
}


Image gradientX(const Image &im, bool clamp) {
  Filter sobelX(3, 3);
  sobelX(0, 0) = -1.0;
  sobelX(1, 0) = 0.0;
  sobelX(2, 0) = 1.0;
  sobelX(0, 1) = -2.0;
  sobelX(1, 1) = 0.0;
  sobelX(2, 1) = 2.0;
  sobelX(0, 2) = -1.0;
  sobelX(1, 2) = 0.0;
  sobelX(2, 2) = 1.0;

  Image imSobelX = sobelX.convolve(im, clamp);
  return imSobelX;
}

Image gradientY(const Image &im, bool clamp) {

  // sobel filtering in y direction
  Filter sobelY(3, 3);
  sobelY(0, 0) = -1.0;
  sobelY(1, 0) = -2.0;
  sobelY(2, 0) = -1.0;
  sobelY(0, 1) = 0.0;
  sobelY(1, 1) = 0.0;
  sobelY(2, 1) = 0.0;
  sobelY(0, 2) = 1.0;
  sobelY(1, 2) = 2.0;
  sobelY(2, 2) = 1.0;

  Image imSobelY = sobelY.convolve(im, clamp);
  return imSobelY;
}

Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
  
  // Compute xx/xy/yy Tensor of an image. (stored in that order)
  //1: compute image dimensions
  int width = im.width();
  int height = im.height();
  int channels = im.channels();

  //2: comput luminance
  vector<Image> lum_chrom = lumiChromi(im);
  Image im_luminance = lum_chrom[0];
  //Image im_chrominance = lum_chrom[1];

  //3: blur using guassian
  //Image gBlur = gaussianBlur_separable(im_luminance,sigmaG);

  //4: compute luminance gradient
  Image Ix = gradientX(im_luminance, true);
  Image Iy = gradientY(im_luminance, true);

  //5: create tensor image
  Image tensor_im(width, height, 3);
  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      float Ix_val = Ix(x, y);
      float Iy_val = Iy(x, y);

      tensor_im(x, y, 0) = Ix_val * Ix_val;
      tensor_im(x, y, 1) = Ix_val * Iy_val;
      tensor_im(x, y, 2) = Iy_val * Iy_val;
    }
  }

  //6: compute sigma tensor
  float sigmaT = sigmaG * factorSigma;
  tensor_im = gaussianBlur_separable(tensor_im, sigmaT);

  //7: return tensor image
  return tensor_im;
}


// Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
//  	// Compute xx/xy/yy Tensor of an image. (stored in that order)

//  	return Image(1,1,1);
// }


Image testAngle(const Image &im, float sigmaG, float factor) {
	// Extracts orientation of features in im. Angles should be mapped
	// to [0,1]
	
	//compute structure tensor 
	Image tensor_im = computeTensor(im, sigmaG, factor);

	//initialise angle map
	Image angleMap(im.width(), im.height(), 1);

	//calculate eigenvalues and eigenvectors for each pixel
	for (int y = 0; y < im.height(); ++y) {
		for (int x = 0; x < im.width(); ++x) {
			//extract tensor components
			float jxx = tensor_im(x, y, 0);
			float jxy = tensor_im(x, y, 1);
			float jyy = tensor_im(x, y, 2);

			//construct structure tensor matrix
			Eigen::Matrix2f tensor;
			tensor(0, 0) = jxx;
			tensor(0, 1) = jxy;
			tensor(1, 0) = jxy;
			tensor(1, 1) = jyy;

			//comput eigenvalues and vectors
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigensolver(tensor);
			if (eigensolver.info() != Eigen::Success) {
				continue;
			}

			//get eigenvalues and eigenvectors
            Eigen::Vector2f eigenvalues = eigensolver.eigenvalues();
            Eigen::Matrix2f eigenvectors = eigensolver.eigenvectors();

            //find the smallest eigenvector
            Eigen::Vector2f smallestEigenvector;
            if (eigenvalues(0) < eigenvalues(1)) {
                smallestEigenvector = eigenvectors.col(0);
            } else {
                smallestEigenvector = eigenvectors.col(1);
            }

            //calculate the angle of the smallest eigenvector with respect to the horizontal axis
            float angle = atan2(smallestEigenvector(1), smallestEigenvector(0));

            // Wrap angle to [0, 2π]
            if (angle < 0) {
                angle += 2 * M_PI;
            }
			

            //map angle to [0, 1] 
            angleMap(x, y, 0) = angle / (2 * M_PI);
			//angleMap(x, y, 0) = 1 - (angle / (2 * M_PI));
		}
	}
    return angleMap;

}

vector<Image> rotateBrushes(const Image &im, int nAngles) {
	// helper function
	// Returns list of nAngles images rotated by 1*2pi/nAngles
	
	vector<Image> rotatedBrushes;
    float angleStep = 2 * M_PI / nAngles;

    for (int i = 0; i < nAngles; ++i) {
        float angle = i * angleStep;
        rotatedBrushes.push_back(rotate(im, angle));  
    }

    return rotatedBrushes;
	
}

void singleScaleOrientedPaint(const Image &im, const Image &importance,
		Image &out, const Image &tensor, const Image &texture,int size, int N, 
		float noise, int nAngles) {
	// Similar to singleScalePaintImportant but brush strokes are oriented
	// according to tensor
	
	//Rotate brushes
    //vector<Image> rotatedBrushes = rotateBrushes(texture, nAngles);

    //Calculate scaling factor for the brush size
    float factor = (texture.width() > texture.height())
                    ? static_cast<float>(size) / texture.width()
                    : static_cast<float>(size) / texture.height();
    
	//scale the texture 
	Image scaledBrush = scaleLin(texture, factor);

	//Rotate scaled brushes
    vector<Image> rotatedBrushes = rotateBrushes(scaledBrush, nAngles);

    int width = out.width();
    int height = out.height();
    int strokesApplied = 0;

    srand(static_cast<unsigned int>(time(0)));

    //apply brush strokes based on orientation and importance sampling
    for (int i = 0; i < N; ++i) {
        int x = rand() % width;
        int y = rand() % height;

        //importance sampling rejection
        float importanceValue = importance(x, y, 0);
        if (static_cast<float>(rand()) / RAND_MAX > importanceValue) {
            continue;  //reject sample
        }

        //extract orientation from tensor at (x, y)
        float jxx = tensor(x, y, 0);
        float jxy = tensor(x, y, 1);
        float jyy = tensor(x, y, 2);

        //construct tensor matrix and find smallest eigenvector angle
        Eigen::Matrix2f tensorMatrix;
        tensorMatrix << jxx, jxy, jxy, jyy;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigensolver(tensorMatrix);
        Eigen::Vector2f smallestEigenvector = eigensolver.eigenvectors().col(0);
        float angle = atan2(smallestEigenvector(1), smallestEigenvector(0));

		//Ensure angle is in [0, 2π]
        if (angle < 0) angle += 2 * M_PI;  

        //map angle to brush rotation index
        int brushIndex = static_cast<int>(round(angle * nAngles / (2 * M_PI))) % nAngles;
        Image &orientedBrush = rotatedBrushes[brushIndex];

        //retrieve color and apply noise factor
        vector<float> color(im.channels());
        float noiseFactor = 1 - noise / 2 + noise * ((float)rand() / RAND_MAX);
        for (int c = 0; c < im.channels(); ++c) {
            color[c] = im(x, y, c) * noiseFactor;
        }

        //apply the oriented brush
        brush(out, x, y, color, orientedBrush);

        strokesApplied++;
        if (strokesApplied >= N) {
            break;
        }
    }
	return;
}

void orientedPaint(const Image &im, Image &out, const Image &texture, int N, int size, float noise) {
	// Similar to painterly() but strokes are oriented along the directions of maximal structure
	
	// Initialize the output image to be a blank canvas
    out = Image(im.width(), im.height(), im.channels());

    //compute the structure tensor for the image
    float sigmaG = 1.0f;
    float factorSigma = 4.0f;
    Image tensor = computeTensor(im, sigmaG, factorSigma);

    //create a constant importance map for the first (coarse) pass
    Image constantImportance(im.width(), im.height(), 1);
    constantImportance.set_color(1.0f);  // Uniform importance for the entire image in the first pass

    //apply the coarse pass using oriented brush strokes
    int coarseSize = size;
    singleScaleOrientedPaint(im, constantImportance, out, tensor, texture, coarseSize, N, noise, 36);
	
    //save intermediate output after the coarse pass
    out.write("./Output/oriented_paint_coarse_pass.png");

    //generate a sharpness map for the second (fine) pass
    float sigma = 1.0f;
    Image sharpness = sharpnessMap(im, sigma);
	//debugging
	sharpness.write("./Output/sharpness_map.png");

    //apply the fine pass using oriented brush strokes, refining details in high sharpness areas
    int fineSize = size / 2;
    int refinedN = N / 2;
    singleScaleOrientedPaint(im, sharpness, out, tensor, texture, fineSize, refinedN, noise, 36);

    //save final output after the fine pass
    out.write("./Output/oriented_paint_final.png");
	

}






//function to predict strokes with parameters based on target image colors, positions, and size
vector<vector<float>> predictStrokes(const Image &target, int numStrokes, float noiseLevel, int size) {
    vector<vector<float>> strokes;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < numStrokes; ++i) {
        vector<float> strokeParams(8);

        //randomly select a position in the target image
        int x = static_cast<int>(dis(gen) * target.width());
        int y = static_cast<int>(dis(gen) * target.height());

        //set position
        strokeParams[0] = x;
        strokeParams[1] = y;

        //set width and height as `size` parameter 
        strokeParams[2] = size;   
        strokeParams[3] = size;    
        strokeParams[4] = dis(gen) * M_PI * 2;  

        //color from target image at (x, y) with added noise
        for (int c = 0; c < 3; ++c) {
            float colorValue = target(x, y, c);
            strokeParams[5 + c] = colorValue * (1.0f - noiseLevel / 2) + noiseLevel * dis(gen);
        }
        
        strokes.push_back(strokeParams);
    }
    return strokes;
}

//helper function to draw a single rectangular stroke directly on the canvas
void drawStroke(Image &canvas, int x, int y, int width, int height, float angle, const vector<float> &color) {
    //apply rotation and draw a filled rectangle representing the stroke
    Image stroke(width, height, 3);
    stroke.set_color(color[0], color[1], color[2]);

    //rotate and place stroke on the canvas
    Image rotatedStroke = rotate(stroke, angle);
    int halfW = rotatedStroke.width() / 2;
    int halfH = rotatedStroke.height() / 2;

    //place rotated stroke on canvas
    for (int j = 0; j < rotatedStroke.height(); ++j) {
        for (int i = 0; i < rotatedStroke.width(); ++i) {
            int canvasX = x - halfW + i;
            int canvasY = y - halfH + j;

            if (canvasX >= 0 && canvasX < canvas.width() && canvasY >= 0 && canvasY < canvas.height()) {
                for (int c = 0; c < canvas.channels(); ++c) {
                    canvas(canvasX, canvasY, c) = rotatedStroke(i, j, c);
                }
            }
        }
    }
}


void paintTransformer(const Image &target, Image &canvas, int numStrokes, float noiseLevel, int numScales) {
    canvas.set_color(1.0f,1.0f,1.0f);  

    for (int scale = 0; scale < numScales; ++scale) {
        int currentSize = 30 / (scale + 1);  

        //generate decreasing strokes
        vector<vector<float>> strokes = predictStrokes(target, numStrokes, noiseLevel, currentSize);

        //render each stroke directly onto the canvas as a shape
        for (const auto &strokeParams : strokes) {
            int x = static_cast<int>(strokeParams[0]);
            int y = static_cast<int>(strokeParams[1]);
            float angle = strokeParams[4];
            vector<float> color = {strokeParams[5], strokeParams[6], strokeParams[7]};

            //draw the stroke as a filled rotated rectangle
            drawStroke(canvas, x, y, currentSize, currentSize, angle, color);
        }


    }

    //save the final output image
    canvas.write("./Output/paint_transformer.png");
}

