/*=========================================================================

  idpMutualInfoVersorRigidVolumeRegistration.txx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef _idpMutualInfoVersorRigidVolumeRegistration_txx
#define _idpMutualInfoVersorRigidVolumeRegistration_txx

#include "idpMutualInfoVersorRigidVolumeRegistration.h"
#include "itkCenteredVersorTransformInitializer.h"

namespace itk {
    namespace idp{

/*
    * Constructor
    */
template < typename TFixedImage, typename TMovingImage >
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::MutualInfoVersorRigidVolumeRegistration() {

    m_FixedInputImage = 0;
    m_MovingInputImage = 0;
    
    // setup default values
    this->SetNumberOfLevels(1);
    m_NumberOfIterations = 100;
    m_MaximumStepLength = 60;
    m_MinimumStepLength = 0.1;
    m_MaximumStepLengthRate = 0.3;
    m_MinimumStepLengthRate = 0.15;
    m_AngleScale = 1.0;
    m_TranslationScale = 1e7;
    m_NumberOfHistogramBins = 64;
    m_NumberOfSpatialSamples = 250000;
    m_RelaxationFactor = 0.8;
    m_RelaxationFactorRate = 1.1;
    m_Verbose = false;

    // setup the default components
    m_OutputTransform = VersorTransformType::New();
    
    typename VersorTransformType::Pointer transform = VersorTransformType::New();
    typename MutualInfoMetricType::Pointer metric = MutualInfoMetricType::New();
    typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    VersorRigid3DTransformOptimizer::Pointer optimizer = VersorRigid3DTransformOptimizer::New();
    m_Coordinator = CoordinatorType::New();

    this->SetTransform(transform);
    this->SetMetric(metric);
    this->SetInterpolator(interpolator);
    this->SetOptimizer(optimizer);

    // set up default values for optimizer
    VersorRigid3DTransformOptimizer::ScalesType optimizerScales(transform->GetNumberOfParameters());
    optimizerScales.Fill(1.0);
    optimizer->SetScales(optimizerScales);
    optimizer->SetMaximize(false); 

    // Setup an optimizer observer
    typedef SimpleMemberCommand<Self> CommandType;
    typename CommandType::Pointer command = CommandType::New();
    command->SetCallbackFunction(this, &Self::OptimizerIterationHandler);
    m_ObserverTag = this->GetGradientDescentOptimizer()->AddObserver(IterationEvent(), command);
}

template < typename TFixedImage, typename TMovingImage >
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::~MutualInfoVersorRigidVolumeRegistration() {
    this->GetGradientDescentOptimizer()->RemoveObserver(m_ObserverTag);
}

template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const {
    this->Superclass::PrintSelf(os, indent);
}

/*
    * Do setup before registration
    */
template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::BeforeRegistration() {
    
    // initialize the coordinator
    m_Coordinator->SetFixedImage( this->m_FixedInputImage );
    m_Coordinator->SetMovingImage( this->m_MovingInputImage );
    m_Coordinator->SetNumberOfLevels( this->GetNumberOfLevels() );
    m_Coordinator->SetFixedImageStartingFactors(8);
    m_Coordinator->SetMovingImageStartingFactors(8);  
    m_Coordinator->Initialize();           
    
    // initialize the transformation
    typedef itk::CenteredVersorTransformInitializer< FixedImageType, MovingImageType>  TransformInitializerType;	
    typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    typename VersorTransformType::Pointer InternalTransform = VersorTransformType::New();
    initializer->SetTransform( InternalTransform );
    initializer->SetFixedImage( m_FixedInputImage );
    initializer->SetMovingImage( m_MovingInputImage );
    initializer->MomentsOn();
    initializer->InitializeTransform();
    this->SetInitialTransformParameters( InternalTransform->GetParameters() );

    // Setup optimizer
    VersorRigid3DTransformOptimizer::Pointer optimizer = this->GetGradientDescentOptimizer();
    optimizer->SetNumberOfIterations(m_NumberOfIterations);
    
    //** set the optimizer 
    typedef VersorRigid3DTransformOptimizer::ScalesType OptimizerScalesType;
    OptimizerScalesType optimizerScales( this->GetTransform()->GetNumberOfParameters() );
    optimizerScales[0] = 1.0 / m_AngleScale; 
    optimizerScales[1] = 1.0 / m_AngleScale; 
    optimizerScales[2] = 1.0 / m_AngleScale; 
    optimizerScales[3] = 1.0 / m_TranslationScale; 
    optimizerScales[4] = 1.0 / m_TranslationScale; 
    optimizerScales[5] = 1.0 / m_TranslationScale; 	
    optimizer->SetScales( optimizerScales );

    optimizer->SetMaximumStepLength( m_MaximumStepLength );  
    optimizer->SetMinimumStepLength( m_MinimumStepLength );

    //*** setting the MI metric
    typename MutualInfoMetricType::Pointer metric = this->GetMutualInfoMetric();
    metric->SetNumberOfHistogramBins( m_NumberOfHistogramBins );
    metric->SetNumberOfSpatialSamples( m_NumberOfSpatialSamples ); 
    metric->ReinitializeSeed( 76926294 );

    optimizer->SetNumberOfIterations( m_NumberOfIterations ); 
    optimizer->SetRelaxationFactor( m_RelaxationFactor );
}

/*
    * Do setup before another iteratio of levels
    */
template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::BeforeIteration() 
{    
    if (m_Verbose) {
        std::cout << "Level: " << this->GetCurrentLevel() << std::endl;
    }

    // Set fixed and moving images for this level
    this->SetFixedImage(m_Coordinator->GetFixedImage(this->GetCurrentLevel()));
    
    typedef typename TFixedImage::RegionType RegionType;
    RegionType region = this->GetFixedImage()->GetBufferedRegion();
    this->SetFixedImageRegion(region);

    this->SetMovingImage(m_Coordinator->GetMovingImage(this->GetCurrentLevel()));

    // Setup optimizer
    if (this->GetCurrentLevel()) {
        typename OptimizerType::Pointer optimizer = this->GetGradientDescentOptimizer();

        // reduce the stepsize in higher level
        optimizer->SetMaximumStepLength( m_MaximumStepLength*m_MaximumStepLengthRate );
        optimizer->SetMinimumStepLength( m_MinimumStepLength*m_MinimumStepLengthRate );
        optimizer->SetRelaxationFactor( m_RelaxationFactor*m_RelaxationFactorRate );
    }
}

/*
    * Do setup after registration
    */
template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::AfterIteration() {

}

/*
    * Do setup after registration
    */
template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::AfterRegistration() {

    // set the output transform
    if ( this->GetLastTransformParameters().Size() != m_OutputTransform->GetNumberOfParameters() )
    {
    itkExceptionMacro(<<"Size mismatch between parameter and transform"); 
    }
    else
    {
    m_OutputTransform->SetParameters( this->GetLastTransformParameters() );
    }

}

/*
    * Print out info each optimization iteration
    */
template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::OptimizerIterationHandler() {
    if (m_Verbose) {
        OptimizerType::Pointer optimizer = this->GetGradientDescentOptimizer();

        std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << optimizer->GetValue() << "   ";
        std::cout << optimizer->GetCurrentPosition() << std::endl;
    }

}

/*
    * Return a pointer to the internal metric
    */
template < typename TFixedImage, typename TMovingImage >
typename MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::MutualInfoMetricType *
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::GetMutualInfoMetric() {
    MutualInfoMetricType * ptr = dynamic_cast<MutualInfoMetricType *>( this->GetMetric() );

    if (!ptr) {
        itkExceptionMacro( << "Metric is not of type NormalizedCorrelationImageToImageMetric");
    }

    return ptr;
}

/*
    * Return a point to the internal transform
    */
template < typename TFixedImage, typename TMovingImage >
typename MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::VersorTransformType *
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::GetVersorTransform() {
    VersorTransformType * ptr =
            dynamic_cast<VersorTransformType *> (this->GetTransform());

    if (!ptr) {
        itkExceptionMacro( << "Metric is not of type Euler3DTransform");
    }

    return ptr;
}

/*
    * Return a point to the internal interpolator
    */
template < typename TFixedImage, typename TMovingImage >
typename MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::LinearInterpolatorType *
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::GetLinearInterpolator() {
    LinearInterpolatorType * ptr =
            dynamic_cast<LinearInterpolatorType *> (this->GetInterpolator());

    if (!ptr) {
        itkExceptionMacro( << "Metric is not of type LinearInterpolateImageFunction");
    }

    return ptr;
}

/*
    * Return a point to the internal optimizer
    */
template < typename TFixedImage, typename TMovingImage >
VersorRigid3DTransformOptimizer *
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>
::GetGradientDescentOptimizer() {
    VersorRigid3DTransformOptimizer * ptr =
            dynamic_cast<VersorRigid3DTransformOptimizer *> (this->GetOptimizer());

    if (!ptr) {
        itkExceptionMacro( << "Metric is not of type GradientDescentOptimizer");
    }

    return ptr;
}

/*
    * read in the registration parameters from a xml file
    */
template < typename TFixedImage, typename TMovingImage >
void
MutualInfoVersorRigidVolumeRegistration<TFixedImage, TMovingImage>::
LoadParametersFromXML( const char * xmlFile )
{
  // open xml file
  TiXmlDocument doc( xmlFile );
  if ( !doc.LoadFile() )
    {
		std::cout<< xmlFile << std::endl;  //debug only
    itkExceptionMacro( << "Could not load xml file " << xmlFile );
    }

  TiXmlNode * node;

  // Locate the root node
  node = doc.FirstChild( "VersorRigid3DRegistration" );
  if ( !node )
    {
    itkExceptionMacro( << "Can not find VersorRigid3DRegistration node in file " << xmlFile );
    }

  TiXmlNode *n;
  TiXmlElement *e;
  n = NULL;

  while( ( n = node->IterateChildren( n ) ) )
    {

    e = n->ToElement();

    if ( e )
      {

      std::string parameter = e->Value();
      if ( parameter.compare( "NumberOfLevels" ) == 0 )
        {
        this->SetNumberOfLevels( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "NumberOfIterations" ) == 0 )
        {
        this->SetNumberOfIterations( atoi( e->GetText() ) );
        }
      else if ( parameter.compare( "MaximumStepLength" ) == 0 )
        {
        this->SetMaximumStepLength( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "MinimumStepLength" ) == 0 )
        {
        this->SetMinimumStepLength( atof( e->GetText() ) );
        }
      else if ( parameter.compare( "MaximumStepLengthRate" ) == 0 )
        {
        this->SetMaximumStepLengthRate( atof( e->GetText() ) );
        }
	  else if ( parameter.compare( "MinimumStepLengthRate" ) == 0 )
        {
        this->SetMinimumStepLengthRate( atof( e->GetText() ) );
        }
          else if ( parameter.compare( "AngleScale" ) == 0 )
        {
        this->SetAngleScale( atof( e->GetText() ) );
        }
	  else if ( parameter.compare( "TranslationScale" ) == 0 )
        {
        this->SetTranslationScale( atof( e->GetText() ) );
        }
	  else if ( parameter.compare( "NumberOfHistogramBins" ) == 0 )
        {
        this->SetNumberOfHistogramBins( atoi( e->GetText() ) );
        }
	  else if ( parameter.compare( "NumberOfSpatialSamples" ) == 0 )
        {
        this->SetNumberOfSpatialSamples( atoi( e->GetText() ) );
        }
	  else if ( parameter.compare( "RelaxationFactor" ) == 0 )
        {
        this->SetRelaxationFactor( atof( e->GetText() ) );
        }
          else if ( parameter.compare( "RelaxationFactorRate" ) == 0 )
        {
        this->SetRelaxationFactorRate( atof( e->GetText() ) );
        }

      } // end if (e)

    } // end while

}

} // end namespace idp

} // end namespace itk


#endif

