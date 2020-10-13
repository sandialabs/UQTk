#include "MUQ/Modeling/WorkPiece.h"

/// An object used to test using user-created inputs to WorkPieces
struct AnObject {
  /// Construct the object
  /**
     @param[in] b The value of the object
   */
  AnObject(double const b) : value(b) {}
  
  /// The object's value
  const double value;
  
  /// A flag that changes the behavior of the object
  bool flag;
};

/// A WorkPiece with no fixed input/output number or type
class UnfixedMod : public muq::Modeling::WorkPiece {
public:

  /// Default constructor
  UnfixedMod() : WorkPiece() {}

  /// Default destructor
  virtual ~UnfixedMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the number of inputs
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
         
    switch( inputs.size() ) {
    case 0 : { // there are no inputs
      const std::string hi = "hello!";
      
      // the output size is 3 
      outputs.resize(3);
      
      outputs[0] = hi;
      outputs[1] = 3.0;
      outputs[2] = 6;
      
      return;
    } case 1: { // there is one input
	
	// the output size is 1
	outputs.resize(1);
	// first input and first ouptut must be a string but the parent does not check this
	outputs[0] = boost::any_cast<std::string>(inputs[0]);
	
	return;
      } case 3: {  // there are three inputs
	  // the first input must be a string (the parent does not check this)
	  const std::string s = boost::any_cast<std::string>(inputs[0]);

	  bool flag = false;
	  if( inputs.size()>1 ) {
	    // the third input must be a shared pointer to AnObject (the parent does not check this)
	    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[2]);
	    flag = obj->flag;
	  }
	  
	  // the output size is 2
	  outputs.resize(2);
	  
	  // the first output is a string
	  outputs[0] = s;
	  
	  if( flag ) { // if the flag in AnObject is true ...
	    // the second input and output must a double (the parent does not check either of these)
	    outputs[1] = boost::any_cast<double>(inputs[1]);
	    
	    return;
	  }
	  
	  // if the flag in AnObject is false
	  outputs[1] = 1.0;
	  
	  return;
	} default: // by default don't do anything      
      return;
    }
  }    
};

/// A WorkPiece with a fixed number of inputs but no fixed input types and no fixed output number or type
class FixedInsMod : public muq::Modeling::WorkPiece {
public:

  /// Constructor 
  /**
     @param[in] numIns The fixed number of inputs
   */
 FixedInsMod(unsigned int numIns) : WorkPiece(numIns) {}

  /// Constructor 
  /**
     @param[in] types The input times
   */
 FixedInsMod(std::vector<std::string> const& types) : WorkPiece(types) {}

  /// Constructor 
  /**
     @param[in] types Some of the input times
   */
 FixedInsMod(std::map<unsigned int, std::string> const& types) : WorkPiece(types) {}
  
  /// Constructor 
  /**
     @param[in] types Some of the input times
     @param[in] numIns The fixed number of inputs
   */
 FixedInsMod(std::map<unsigned int, std::string> const& types, unsigned int const numIns) : WorkPiece(types, numIns) {}

  /// Default destructor
  virtual ~FixedInsMod() {}

private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the type of the third input
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the third input must be a shared pointer to AnObject (the parent may check this)
    auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[2]);

    if( obj->flag ) { // if the third input is a bool
      // there are 2 outputs
      outputs.resize(2);

      // first input and first output must be a string (the parent may check this)
      outputs[0] = boost::any_cast<std::string>(inputs[0]);
      // second input and second output must be an unsigned int (the parent may check this)
      outputs[1] = obj->value*boost::any_cast<double>(inputs[1]);
    } else { // if the third input is not a bool
      // there is 1 output
      outputs.resize(1);

      // first output and second input must be an unsigned int (the parent may check this)
      outputs[0] = boost::any_cast<double>(inputs[1]);
    }
  }    
};

/// A WorkPiece with no fixed input number or type and a fixed number of ouputs but no fixed output type
class FixedOutsMod : public muq::Modeling::WorkPiece {
 public:
  
  /**
     @param[in] numOuts The fixed number of outputs
  */
 FixedOutsMod(unsigned int numOuts) : WorkPiece(numOuts, WorkPiece::Fix::Outputs) {}
  
  /**
     @param[in] types The output times
  */
 FixedOutsMod(std::vector<std::string> const& types) : WorkPiece(types, WorkPiece::Fix::Outputs) {}
  
  /**
     @param[in] types Some of the output times
  */
 FixedOutsMod(std::map<unsigned int, std::string> const& types) : WorkPiece(types, WorkPiece::Fix::Outputs) {}
  
  /**
     @param[in] types Some of the output times
     @param[in] num The number of outputs
  */
 FixedOutsMod(std::map<unsigned int, std::string> const& types, unsigned int const num) : WorkPiece(types, num, WorkPiece::Fix::Outputs, WorkPiece::Fix::Outputs) {}
  
  /// Default destructor
  virtual ~FixedOutsMod() {}
  
 private:

  /// User-defined EvaluateImpl function
  /**
     The behavior changes depending on the number of inputs
   */
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string (the parent may check this)
    const std::string s = inputs.size()>0? boost::any_cast<std::string>(inputs[0]) : "";
    // the second input must be a shared pointer to an object (the parent may check this)
    auto obj = inputs.size()>1? boost::any_cast<std::shared_ptr<AnObject> >(inputs[1]) : nullptr;

    outputs.resize(2);

    // the first output is a string (the parent may check this)
    outputs[0] = s;

    // the second object is a double (the parent may check this)
    outputs[1] = obj? obj->value : 1.0;
  }    
};

/// A WorkPiece with a fixed number of inputs and outputs but no fixed types
class FixedInOutMod : public muq::Modeling::WorkPiece {
public:

  /**
     @param[in] numIns The fixed number of inputs
     @param[in] numOuts The fixed number of outputs
   */
  FixedInOutMod(unsigned int const numIns, unsigned int numOuts) : WorkPiece(numIns, numOuts) {}

  /**
     @param[in] inTypes Some of the input times
     @param[in] numOuts The fixed number of outputs
   */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, unsigned int numOuts) : WorkPiece(inTypes, numOuts, WorkPiece::Fix::Inputs, WorkPiece::Fix::Outputs) {}

  /**
     @param[in] inTypes The input times
     @param[in] numOuts The fixed number of outputs
  */
 FixedInOutMod(std::vector<std::string> const& inTypes, unsigned int numOuts) : WorkPiece(inTypes, numOuts) {}

  /**
     @param[in] numIns The fixed number of inputs
     @param[in] outTypes The output times
  */
 FixedInOutMod(unsigned int const numIns, std::vector<std::string> const& outTypes) : WorkPiece(numIns, outTypes) {}
  
  /**
     @param[in] numIns The fixed number of inputs
     @param[in] outTypes Some of the output times
  */
 FixedInOutMod(unsigned int numIns, std::map<unsigned int, std::string> const& outTypes) : WorkPiece(outTypes, numIns, WorkPiece::Fix::Outputs) {}

  /**
     @param[in] types The input times
     @param[in] inNum The number of inputs
     @param[in] outNum The number of outputs
   */
 FixedInOutMod(std::map<unsigned int, std::string> const& types, unsigned int const inNum, unsigned int const outNum) : WorkPiece(types, inNum, outNum) {}

  /**
     @param[in] inNum The number of inputs
     @param[in] types The output times
     @param[in] outNum The number of outputs
   */
 FixedInOutMod(unsigned int const inNum, std::map<unsigned int, std::string> const& types, unsigned int const outNum) : WorkPiece(inNum, types, outNum) {}

  /**
     @param[in] inTypes The input types
     @param[in] outTypes The output types
   */
 FixedInOutMod(std::vector<std::string> const& inTypes, std::vector<std::string> const& outTypes) : WorkPiece(inTypes, outTypes) {}

  /**
     @param[in] inTypes Some of the input types
     @param[in] outTypes The output types
   */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, std::vector<std::string> const& outTypes) : WorkPiece(inTypes, outTypes) {}

  /**
     @param[in] inTypes Some of the input types
     @param[in] numIns The number of inputs
     @param[in] outTypes The output types
   */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, unsigned int const numIns, std::vector<std::string> const& outTypes) : WorkPiece(inTypes, numIns, outTypes) {}

  /**
     @param[in] inTypes The input types
     @param[in] outTypes Some of the output types
  */
 FixedInOutMod(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes) : WorkPiece(inTypes, outTypes) {}

  /**
     @param[in] inTypes The input types
     @param[in] outTypes Some of the output types
     @param[in] numOut The number of outputs
  */
 FixedInOutMod(std::vector<std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, unsigned int const numOut) : WorkPiece(inTypes, outTypes, numOut) {}

  /**
     @param[in] inTypes Some of the input types
     @param[in] outTypes Some of the output types
  */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes) : WorkPiece(inTypes, outTypes) {}

  /**
     @param[in] inTypes Some of the input types
     @param[in] numIns The number of inputs
     @param[in] outTypes Some of the output types
  */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, unsigned int const numIns, std::map<unsigned int, std::string> const& outTypes) : WorkPiece(inTypes, numIns, outTypes) {}

    /**
     @param[in] inTypes Some of the input types
     @param[in] outTypes Some of the output types
     @param[in] numOuts The number of outputs
  */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, std::map<unsigned int, std::string> const& outTypes, unsigned int const numOuts) : WorkPiece(inTypes, outTypes, numOuts) {}

    /**
     @param[in] inTypes Some of the input types
     @param[in] numIns The number of inputs
     @param[in] outTypes Some of the output types
     @param[in] numOuts The number of outputs
  */
 FixedInOutMod(std::map<unsigned int, std::string> const& inTypes, unsigned int const numIns, std::map<unsigned int, std::string> const& outTypes, unsigned int const numOuts) : WorkPiece(inTypes, numIns, outTypes, numOuts) {}

  /// Default destructor
  virtual ~FixedInOutMod() {}

private:

  /// User-defined EvaluateImpl function
  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    // the first input must be a string (the parent may check this)
    const std::string s = boost::any_cast<std::string>(inputs[0]);

    bool flag = false;
    if( inputs.size()>1 ) {
      // the third input must be a shared pointer to AnObject (the parent may check this)
      auto obj = boost::any_cast<std::shared_ptr<AnObject> >(inputs[2]);
      flag = obj->flag;
    }

    // the output size is 2
    outputs.resize(2);

    // the first output is a string
    outputs[0] = s;

    if( flag ) { // if the flag in AnObject is true ...
      // the second input and output must a double (the parent may check either or both of these)
      outputs[1] = boost::any_cast<double>(inputs[1]);

      return;
    }

    // if the flag in AnObject is false
    outputs[1] = 1.0;
  }
};

