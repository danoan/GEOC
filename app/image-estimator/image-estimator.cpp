#include <unistd.h>

#include <DGtal/helpers/StdDefs.h>
#include <opencv2/highgui.hpp>
#include <boost/filesystem.hpp>

#include <DIPaCUS/base/Shapes.h>
#include <DIPaCUS/derivates/Misc.h>
#include <DIPaCUS/components/Transform.h>
#include <geoc/api/api.h>

struct InputData
{
    enum class Estimator{ISC_MDCA,ISC_II,SQC_MDCA,SQC_II,LENGTH_PROJECTION,LENGTH_SIN};
    enum class Mode{AllInFolder,SingleImage,Shape};
    enum class Shape{Triangle,Square,Flower,Ball,Bean};

    InputData()
    {
        inputPath="";
        outputPath="";

        mode=Mode::Shape;
        shape=Shape::Triangle;
        estimator=Estimator::ISC_MDCA;

        gridStep = 1.0;
        radius = 5.0;
        lengthPenalization=0;
    }

    std::string resolve(Mode mode) const
    {
        switch(mode)
        {
            case Mode::Shape: return "shape";
            case Mode::AllInFolder: return "all-in-folder";
            case Mode::SingleImage: return "single-image";
            default: return "Not recognized:";
        }
    }

    std::string resolve(Estimator estimator) const
    {
        switch(estimator)
        {
            case Estimator::ISC_MDCA: return "mdca";
            case Estimator::ISC_II: return "ii";
            case Estimator::LENGTH_PROJECTION: return "length-projection";
            case Estimator::LENGTH_SIN: return "length-sin";
            default: return "Not recognized:";
        }
    }

    std::string resolve(Shape shape) const
    {
        switch(shape)
        {
            case Shape::Triangle: return "triangle";
            case Shape::Square: return "square";
            case Shape::Flower: return "flower";
            case Shape::Ball: return "ball";
            case Shape::Bean: return "bean";
            default: return "Not recognized:";
        }
    }

    std::string inputPath;
    std::string outputPath;

    Mode mode;
    Shape shape;
    Estimator estimator;

    double gridStep;
    double radius;
    double lengthPenalization;
};

void usage(int argc, char* argv[])
{
    std::cerr << "Usage " << argv[0] << ":\n"
    << "[-m] Mode (shape,all-in-folder, single-image)\n"
    << "[-s] Shape (triangle square flower ball bean)\n"
    << "[-f] Image/Folder path\n"
    << "[-e] Estimator (isc-mdca,isc-ii,sqc-mdca,sqc-ii,length-projection,length-sin)\n"
    << "[-h] Grid Step\n"
    << "[-r] II estimation ball radius\n"
    << "[-a] Length penalization (SQC)\n"
    << "[-o] OutputFilepath \n";
}

InputData readInput(int argc, char* argv[])
{
    if(argc<2)
    {
        usage(argc,argv);
        exit(0);
    }

    InputData id;

    int opt;
    while( ( opt=getopt(argc,argv,"m:s:f:e:h:o:r:a:") )!=-1 )
    {
        switch(opt)
        {
            case 'm':
            {
                if(strcmp(optarg,"shape")==0) id.mode = InputData::Mode::Shape;
                else if(strcmp(optarg,"all-in-folder")==0) id.mode = InputData::Mode::AllInFolder;
                else if(strcmp(optarg,"single-image")==0) id.mode = InputData::Mode::SingleImage;
                else throw std::runtime_error("Mode not recognized.");
                break;
            }
            case 's':
            {
                if(strcmp(optarg,"triangle")==0) id.shape = InputData::Shape::Triangle;
                else if(strcmp(optarg,"square")==0) id.shape = InputData::Shape::Square;
                else if(strcmp(optarg,"flower")==0) id.shape = InputData::Shape::Flower;
                else if(strcmp(optarg,"ball")==0) id.shape = InputData::Shape::Ball;
                else if(strcmp(optarg,"bean")==0) id.shape = InputData::Shape::Bean;
                else throw std::runtime_error("Shape not recognized.");
                break;
            }
            case 'f':
            {
                id.inputPath=optarg;
                break;
            }
            case 'e':
            {
                if(strcmp(optarg,"isc-mdca")==0) id.estimator= InputData::Estimator::ISC_MDCA;
                else if(strcmp(optarg,"isc-ii")==0) id.estimator = InputData::Estimator::ISC_II;
                else if(strcmp(optarg,"sqc-mdca")==0) id.estimator = InputData::Estimator::SQC_MDCA;
                else if(strcmp(optarg,"sqc-ii")==0) id.estimator = InputData::Estimator::SQC_II;
                else if(strcmp(optarg,"length-projection")==0) id.estimator = InputData::Estimator::LENGTH_PROJECTION;
                else if(strcmp(optarg,"length-sin")==0) id.estimator = InputData::Estimator::LENGTH_SIN;
                else throw std::runtime_error("Estimator not recognized.");
                break;
            }
            case 'h':
            {
                id.gridStep=std::atof(optarg);
                break;
            }
            case 'r':
            {
                id.radius=std::atof(optarg);
                break;
            }
            case 'a':
            {
                id.lengthPenalization=std::atof(optarg);
                break;
            }
            case 'o':
            {
                id.outputPath = optarg;
                break;
            }
            default:
            {
                usage(argc,argv);
                exit(0);
            }
        }
    }

    return id;
}

typedef DGtal::Z2i::DigitalSet DigitalSet;
typedef DGtal::Z2i::Domain Domain;
typedef DGtal::Z2i::Point Point;

DigitalSet resolveShape(InputData::Shape shape,double gridStep)
{
    switch(shape)
    {
        case InputData::Shape::Triangle: return DIPaCUS::Shapes::triangle(gridStep);
        case InputData::Shape::Square: return DIPaCUS::Shapes::square(gridStep);
        case InputData::Shape::Flower: return DIPaCUS::Shapes::flower(gridStep);
        case InputData::Shape::Ball: return DIPaCUS::Shapes::ball(gridStep);
        case InputData::Shape::Bean: return DIPaCUS::Shapes::bean(gridStep);
        default: throw std::runtime_error("Shape not recognized!");
    }
}

double runEstimation(const DigitalSet& _ds, InputData::Estimator& estimator, double gridStep, double radius,double lengthPenalization)
{
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

    Curve curve;
    std::vector<double> evK;
    std::vector<double> evS;

    DigitalSet ds = DIPaCUS::Misc::cleanSet(_ds);

    DIPaCUS::Misc::computeBoundaryCurve(curve,ds);
    KSpace kspace;
    kspace.init(ds.domain().lowerBound(),ds.domain().upperBound(),true);

    GEOC::Estimator::Standard::IICurvatureExtraData iiData(true,radius);
    double v=0;
    switch(estimator)
    {
        case InputData::Estimator::ISC_MDCA:
        {
            using namespace GEOC::API::GridCurve;
            Curvature::symmetricClosed<Curvature::EstimationAlgorithms::ALG_MDCA>(kspace,curve.begin(),curve.end(),evK,gridStep,&iiData);
            Length::mdssClosed<Length::EstimationAlgorithms::ALG_PROJECTED>(kspace,curve.begin(),curve.end(),evS,gridStep,NULL);

            for( auto i=0;i<evK.size();++i ) v+= pow(evK[i],2)*evS[i];
            break;
        }
        case InputData::Estimator::ISC_II:
        {
            using namespace GEOC::API::GridCurve;
            Curvature::identityOpen<Curvature::EstimationAlgorithms::ALG_II>(kspace,curve.begin(),curve.end(),evK,gridStep,&iiData);
            Length::mdssClosed<Length::EstimationAlgorithms::ALG_PROJECTED>(kspace,curve.begin(),curve.end(),evS,gridStep,NULL);

            for( auto i=0;i<evK.size();++i ) v+= pow(evK[i],2)*evS[i];
            break;
        }
        case InputData::Estimator::SQC_MDCA:
        {
            using namespace GEOC::API::GridCurve;
            Curvature::symmetricClosed<Curvature::EstimationAlgorithms::ALG_MDCA>(kspace,curve.begin(),curve.end(),evK,gridStep,&iiData);
            Length::mdssClosed<Length::EstimationAlgorithms::ALG_PROJECTED>(kspace,curve.begin(),curve.end(),evS,gridStep,NULL);

            for( auto i=0;i<evK.size();++i ) v+= gridStep*pow(evK[i],2) + lengthPenalization*evS[i];
            break;
        }
        case InputData::Estimator::SQC_II:
        {
            using namespace GEOC::API::GridCurve;
            Curvature::identityOpen<Curvature::EstimationAlgorithms::ALG_II>(kspace,curve.begin(),curve.end(),evK,gridStep,&iiData);
            Length::mdssClosed<Length::EstimationAlgorithms::ALG_PROJECTED>(kspace,curve.begin(),curve.end(),evS,gridStep,NULL);

            for( auto i=0;i<evK.size();++i ) v+= gridStep*pow(evK[i],2) + lengthPenalization*evS[i];
            break;
        }
        case InputData::Estimator::LENGTH_PROJECTION:
        {
            using namespace GEOC::API::GridCurve::Length;
            mdssClosed<EstimationAlgorithms::ALG_PROJECTED>(kspace,curve.begin(),curve.end(),evS,gridStep,NULL);
            for( double x:evS ) v+=x;
            break;
        }
        case InputData::Estimator::LENGTH_SIN:
        {
            using namespace GEOC::API::GridCurve::Length;
            mdssClosed<EstimationAlgorithms::ALG_SINCOS>(kspace,curve.begin(),curve.end(),evS,gridStep,NULL);
            for( double x:evS ) v+=x;
            break;
        }
    }



    return v;
}

DigitalSet digitalSetFromImagePath(const std::string& imagePath)
{
    auto img = cv::imread(imagePath,DIPaCUS::Representation::GRAYSCALE_IMG_TYPE);
    Domain domain( Point(0,0), Point(img.cols-1,img.rows-1) );
    DigitalSet ds(domain);
    DIPaCUS::Representation::CVMatToDigitalSet(ds,img);

    return ds;
}

void writeInputData(std::ostream& ofs, const InputData& id)
{
    ofs << "Input Path: " << id.inputPath << "\n"
    << "Output Path: " << id.outputPath << "\n"
    << "Mode: " << id.resolve(id.mode) << "\n"
    << "Shape: " << id.resolve(id.shape) << "\n"
    << "Estimator: " << id.resolve(id.estimator) << "\n"
    << "Grid Step: " << id.gridStep << "\n"
    << "Radius: " << id.radius << "\n\n";
}

void writeEstimationValue(std::ostream& ofs, const std::string& name, double value)
{
    ofs << name << ": " << value << "\n";
}

int main(int argc, char* argv[])
{
    InputData id = readInput(argc,argv);

    std::ostream* ofs;
    if(id.outputPath=="") ofs = &std::cout;
    else
    {
        boost::filesystem::path p(id.outputPath);
        boost::filesystem::create_directories(p.remove_filename());

        ofs = new std::ofstream(id.outputPath);
        writeInputData(*ofs,id);
    }



    double v;
    switch(id.mode)
    {
        case InputData::Mode::Shape:
        {
            v=runEstimation(resolveShape(id.shape,id.gridStep),id.estimator,id.gridStep,id.radius,id.lengthPenalization);
            writeEstimationValue(*ofs,"Shape",v);
            break;
        }
        case InputData::Mode::SingleImage:
        {
            v=runEstimation(digitalSetFromImagePath(id.inputPath),id.estimator,id.gridStep,id.radius,id.lengthPenalization);
            writeEstimationValue(*ofs,"Single-Image",v);
            break;
        }
        case InputData::Mode::AllInFolder:
        {
            using namespace boost::filesystem;
            path p(id.inputPath);
            directory_iterator di(p);
            while(di!=directory_iterator())
            {
                if(is_regular_file(*di))
                {
                    path curr_path = di->path();
                    v=runEstimation(digitalSetFromImagePath(curr_path.string()),id.estimator,id.gridStep,id.radius,id.lengthPenalization);
                    writeEstimationValue(*ofs,curr_path.filename().string(),v);
                }
                ++di;
            }
            break;
        }
    }

    if(id.outputPath!="")
    {
        ofs->flush();
        delete ofs;
    }

    return 0;
}