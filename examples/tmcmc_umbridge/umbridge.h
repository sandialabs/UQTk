#ifndef UMBRIDGE
#define UMBRIDGE

// Increase timeout to allow for long-running models.
// This should be (to be on the safe side) significantly greater than the maximum time your model may take
#define CPPHTTPLIB_READ_TIMEOUT_SECOND 60*60

#include <string>
#include <vector>


#include "json.hpp"
#include "httplib.h"

using json = nlohmann::json;

namespace umbridge {

  class Model {
  public:
    Model(std::string name) : name(name) {}

    virtual std::vector<std::size_t> GetInputSizes(const json& config_json = json::parse("{}")) const = 0;
    virtual std::vector<std::size_t> GetOutputSizes(const json& config_json = json::parse("{}")) const = 0;

    virtual std::vector<std::vector<double>> Evaluate(const std::vector<std::vector<double>>& inputs,
                          json config_json = json::parse("{}")) {
      (void)inputs; (void)config_json; // Avoid unused argument warnings
      throw std::runtime_error("Evaluate was called, but not implemented by model!");
    }

    virtual std::vector<double> Gradient(unsigned int outWrt,
                          unsigned int inWrt,
                          const std::vector<std::vector<double>>& inputs,
                          const std::vector<double>& sens,
                          json config_json = json::parse("{}")) {
      (void)outWrt; (void)inWrt; (void)inputs; (void)sens; (void)config_json; // Avoid unused argument warnings
      throw std::runtime_error("Gradient was called, but not implemented by model!");
    }

    virtual std::vector<double> ApplyJacobian(unsigned int outWrt,
                              unsigned int inWrt,
                              const std::vector<std::vector<double>>& inputs,
                              const std::vector<double>& vec,
                              json config_json = json::parse("{}")) {
      (void)outWrt; (void)inWrt; (void)inputs; (void)vec; (void)config_json; // Avoid unused argument warnings
      throw std::runtime_error("ApplyJacobian was called, but not implemented by model!");
    }

    virtual std::vector<double> ApplyHessian(unsigned int outWrt,
                              unsigned int inWrt1,
                              unsigned int inWrt2,
                              const std::vector<std::vector<double>>& inputs,
                              const std::vector<double>& sens,
                              const std::vector<double>& vec,
                              json config_json = json::parse("{}")) {
      (void)outWrt; (void)inWrt1; (void)inWrt2; (void)inputs; (void)sens; (void)vec; (void)config_json; // Avoid unused argument warnings
      throw std::runtime_error("ApplyHessian was called, but not implemented by model!");
    }

    virtual bool SupportsEvaluate() {return false;}
    virtual bool SupportsGradient() {return false;}
    virtual bool SupportsApplyJacobian() {return false;}
    virtual bool SupportsApplyHessian() {return false;}

    std::string GetName() const {return name;}

  protected:
    std::string name;
  };

  std::vector<std::string> SupportedModels(std::string host, httplib::Headers headers = httplib::Headers()) {
    httplib::Client cli(host.c_str());
    if (auto res = cli.Get("/Info", headers)) {
      json response = json::parse(res->body);

      if (response.value<double>("protocolVersion",0) != 1.0)
        throw std::runtime_error("Model protocol version not supported!");

      return response["models"];

    } else {
      throw std::runtime_error("GET Info failed with error type '" + to_string(res.error()) + "'");
    }
  }

  // Client-side Model connecting to a server for the actual evaluations etc.
  class HTTPModel : public Model {
  public:

    HTTPModel(std::string host, std::string name, httplib::Headers headers = httplib::Headers())
    : Model(name), cli(host.c_str()), headers(headers)
    {
      // Check if requested model is available on server
      std::vector<std::string> models = SupportedModels(host, headers);
      if (std::find(models.begin(), models.end(), name) == models.end()) {
        std::string model_names = "";
        for (auto& m : models) {
          model_names += "'" + m + "' ";
        }
        throw std::runtime_error("Model " + name + " not found on server! Available models: " + model_names + ".");
      }

      json request_body;
      request_body["name"] = name;

      if (auto res = cli.Post("/ModelInfo", headers, request_body.dump(), "application/json")) {
        json response = json::parse(res->body);

        json supported_features = response.at("support");
        supportsEvaluate = supported_features.value("Evaluate", false);
        supportsGradient = supported_features.value("Gradient", false);
        supportsApplyJacobian = supported_features.value("ApplyJacobian", false);
        supportsApplyHessian = supported_features.value("ApplyHessian", false);
      } else {
        throw std::runtime_error("POST ModelInfo failed with error type '" + to_string(res.error()) + "'");
      }
    }

    std::vector<std::size_t> GetInputSizes(const json& config_json = json::parse("{}")) const override {

      json request_body;
      request_body["name"] = name;
      if (!config_json.empty())
        request_body["config"] = config_json;

      if (auto res = cli.Post("/InputSizes", headers, request_body.dump(), "application/json")) {
        json response_body = parse_result_with_error_handling(res);
        std::vector<std::size_t> outputvec = response_body["inputSizes"].get<std::vector<std::size_t>>();
        return outputvec;
      } else {
        throw std::runtime_error("POST InputSizes failed with error type '" + to_string(res.error()) + "'");
        return std::vector<std::size_t>(0);
      }
    }

    std::vector<std::size_t> GetOutputSizes(const json& config_json = json::parse("{}")) const override {

      json request_body;
      request_body["name"] = name;
      if (!config_json.empty())
        request_body["config"] = config_json;

      if (auto res = cli.Post("/OutputSizes", headers, request_body.dump(), "application/json")) {
        json response_body = parse_result_with_error_handling(res);
        std::vector<std::size_t> outputvec = response_body["outputSizes"].get<std::vector<std::size_t>>();
        return outputvec;
      } else {
        throw std::runtime_error("POST OutputSizes failed with error type '" + to_string(res.error()) + "'");
        return std::vector<std::size_t>(0);
      }
    }

    std::vector<std::vector<double>> Evaluate(const std::vector<std::vector<double>>& inputs, json config_json = json::parse("{}")) override {

      json request_body;
      request_body["name"] = name;

      for (std::size_t i = 0; i < inputs.size(); i++) {
        request_body["input"][i] = inputs[i];
      }
      request_body["config"] = config_json;

      if (auto res = cli.Post("/Evaluate", headers, request_body.dump(), "application/json")) {
        json response_body = parse_result_with_error_handling(res);

        std::vector<std::vector<double>> outputs(response_body["output"].size());
        for (std::size_t i = 0; i < response_body["output"].size(); i++) {
          outputs[i] = response_body["output"][i].get<std::vector<double>>();
        }
        return outputs;
      } else {
        throw std::runtime_error("POST Evaluate failed with error type '" + to_string(res.error()) + "'");
      }
    }

    std::vector<double> Gradient(unsigned int outWrt,
                  unsigned int inWrt,
                  const std::vector<std::vector<double>>& inputs,
                  const std::vector<double>& sens,
                  json config_json = json::parse("{}")) override
    {

      json request_body;
      request_body["name"] = name;
      request_body["outWrt"] = outWrt;
      request_body["inWrt"] = inWrt;
      for (std::size_t i = 0; i < inputs.size(); i++) {
        request_body["input"][i] = inputs[i];
      }
      request_body["sens"] = sens;
      request_body["config"] = config_json;

      if (auto res = cli.Post("/Gradient", headers, request_body.dump(), "application/json")) {
        json response_body = parse_result_with_error_handling(res);

        return response_body["output"].get<std::vector<double>>();
      } else {
        throw std::runtime_error("POST Gradient failed with error type '" + to_string(res.error()) + "'");
      }
    }

    std::vector<double> ApplyJacobian(unsigned int outWrt,
                              unsigned int inWrt,
                              const std::vector<std::vector<double>>& inputs,
                              const std::vector<double>& vec,
                              json config_json = json::parse("{}")) override {

      json request_body;
      request_body["name"] = name;
      request_body["outWrt"] = outWrt;
      request_body["inWrt"] = inWrt;
      for (std::size_t i = 0; i < inputs.size(); i++) {
        request_body["input"][i] = inputs[i];
      }
      request_body["vec"] = vec;
      request_body["config"] = config_json;

      if (auto res = cli.Post("/ApplyJacobian", headers, request_body.dump(), "application/json")) {
        json response_body = parse_result_with_error_handling(res);

        return response_body["output"].get<std::vector<double>>();
      } else {
        throw std::runtime_error("POST ApplyJacobian failed with error type '" + to_string(res.error()) + "'");
      }
    }

    std::vector<double> ApplyHessian(unsigned int outWrt,
                      unsigned int inWrt1,
                      unsigned int inWrt2,
                      const std::vector<std::vector<double>>& inputs,
                      const std::vector<double>& sens,
                      const std::vector<double>& vec,
                      json config_json = json::parse("{}")) override {

      json request_body;
      request_body["name"] = name;
      request_body["outWrt"] = outWrt;
      request_body["inWrt1"] = inWrt1;
      request_body["inWrt2"] = inWrt2;
      for (std::size_t i = 0; i < inputs.size(); i++) {
        request_body["input"][i] = inputs[i];
      }
      request_body["sens"] = sens;
      request_body["vec"] = vec;
      request_body["config"] = config_json;

      if (auto res = cli.Post("/ApplyHessian", headers, request_body.dump(), "application/json")) {
        json response_body = parse_result_with_error_handling(res);

        return response_body["output"].get<std::vector<double>>();
      } else {
        throw std::runtime_error("POST ApplyHessian failed with error type '" + to_string(res.error()) + "'");
      }
    }

    bool SupportsEvaluate() override {
      return supportsEvaluate;
    }
    bool SupportsGradient() override {
      return supportsGradient;
    }
    bool SupportsApplyJacobian() override {
      return supportsApplyJacobian;
    }
    bool SupportsApplyHessian() override {
      return supportsApplyHessian;
    }

  private:

    mutable httplib::Client cli;
    httplib::Headers headers;

    bool supportsEvaluate = false;
    bool supportsGradient = false;
    bool supportsApplyJacobian = false;
    bool supportsApplyHessian = false;

    json parse_result_with_error_handling(const httplib::Result& res) const {
      json response_body;
      try {
        response_body = json::parse(res->body);
      } catch (json::parse_error& e) {
        throw std::runtime_error("Response JSON could not be parsed. Response body: '" + res->body + "'");
      }
      if (response_body.find("error") != response_body.end()) {
        throw std::runtime_error("Model server returned error of type " + response_body["error"]["type"].get<std::string>() + ", message: " + response_body["error"]["message"].get<std::string>());
      }
      return response_body;
    }

  };

  // Check if inputs dimensions match model's expected input size and return error in httplib response
  bool check_input_sizes(const std::vector<std::vector<double>>& inputs, const json& config_json, const Model& model, httplib::Response& res) {
    if (inputs.size() != model.GetInputSizes(config_json).size()) {
      json response_body;
      response_body["error"]["type"] = "InvalidInput";
      response_body["error"]["message"] = "Number of inputs does not match number of model inputs. Expected " + std::to_string(model.GetInputSizes(config_json).size()) + " but got " + std::to_string(inputs.size());
      res.set_content(response_body.dump(), "application/json");
      res.status = 400;
      return false;
    }
    for (std::size_t i = 0; i < inputs.size(); i++) {
      if (inputs[i].size() != model.GetInputSizes(config_json)[i]) {
        json response_body;
        response_body["error"]["type"] = "InvalidInput";
        response_body["error"]["message"] = "Input size mismatch! In input " + std::to_string(i) + " model expected size " + std::to_string(model.GetInputSizes(config_json)[i]) + " but got " + std::to_string(inputs[i].size());
        res.set_content(response_body.dump(), "application/json");
        res.status = 400;
        return false;
      }
    }
    return true;
  }

  // Check if sensitivity vector's dimension matches correct model output size and return error in httplib response
  bool check_sensitivity_size(const std::vector<double>& sens, int outWrt, const json& config_json, const Model& model, httplib::Response& res) {
    if (sens.size() != model.GetOutputSizes(config_json)[outWrt]) {
      json response_body;
      response_body["error"]["type"] = "InvalidInput";
      response_body["error"]["message"] = "Sensitivity vector size mismatch! Expected " + std::to_string(model.GetOutputSizes(config_json)[outWrt]) + " but got " + std::to_string(sens.size());
      res.set_content(response_body.dump(), "application/json");
      res.status = 400;
      return false;
    }
    return true;
  }

  // Check if vector's dimension matches correct model output size and return error in httplib response
  bool check_vector_size(const std::vector<double>& vec, int inWrt, const json& config_json, const Model& model, httplib::Response& res) {
    if (vec.size() != model.GetInputSizes(config_json)[inWrt]) {
      json response_body;
      response_body["error"]["type"] = "InvalidInput";
      response_body["error"]["message"] = "Vector size mismatch! Expected " + std::to_string(model.GetInputSizes(config_json)[inWrt]) + " but got " + std::to_string(vec.size());
      res.set_content(response_body.dump(), "application/json");
      res.status = 400;
      return false;
    }
    return true;
  }

  // Check if outputs dimensions match model's expected output size and return error in httplib response
  bool check_output_sizes(const std::vector<std::vector<double>>& outputs, const json& config_json, const Model& model, httplib::Response& res) {
    if (outputs.size() != model.GetOutputSizes(config_json).size()) {
      json response_body;
      response_body["error"]["type"] = "InvalidOutput";
      response_body["error"]["message"] = "Number of outputs declared by model does not match number of outputs returned by model. Model declared " + std::to_string(model.GetOutputSizes(config_json).size()) + " but returned " + std::to_string(outputs.size());
      res.set_content(response_body.dump(), "application/json");
      res.status = 500;
      return false;
    }
    for (std::size_t i = 0; i < outputs.size(); i++) {
      if (outputs[i].size() != model.GetOutputSizes(config_json)[i]) {
        json response_body;
        response_body["error"]["type"] = "InvalidOutput";
        response_body["error"]["message"] = "Output size mismatch! In output " + std::to_string(i) + " model declared size " + std::to_string(model.GetOutputSizes(config_json)[i]) + " but returned " + std::to_string(outputs[i].size());
        res.set_content(response_body.dump(), "application/json");
        res.status = 500;
        return false;
      }
    }
    return true;
  }

  // Check if inWrt is between zero and model's input size inWrt and return error in httplib response
  bool check_input_wrt(int inWrt, const json& config_json, const Model& model, httplib::Response& res) {
    if (inWrt < 0 || inWrt >= (int)model.GetInputSizes(config_json).size()) {
      json response_body;
      response_body["error"]["type"] = "InvalidInput";
      response_body["error"]["message"] = "Input inWrt out of range! Expected between 0 and " + std::to_string(model.GetInputSizes(config_json).size() - 1) + " but got " + std::to_string(inWrt);
      res.set_content(response_body.dump(), "application/json");
      res.status = 400;
      return false;
    }
    return true;
  }

  // Check if outWrt is between zero and model's output size outWrt and return error in httplib response
  bool check_output_wrt(int outWrt, const json& config_json, const Model& model, httplib::Response& res) {
    if (outWrt < 0 || outWrt >= (int)model.GetOutputSizes(config_json).size()) {
      json response_body;
      response_body["error"]["type"] = "InvalidInput";
      response_body["error"]["message"] = "Input outWrt out of range! Expected between 0 and " + std::to_string(model.GetOutputSizes(config_json).size() - 1) + " but got " + std::to_string(outWrt);
      res.set_content(response_body.dump(), "application/json");
      res.status = 400;
      return false;
    }
    return true;
  }

  // Construct response for unsupported feature
  void write_unsupported_feature_response(httplib::Response& res, std::string feature) {
    json response_body;
    response_body["error"]["type"] = "UnsupportedFeature";
    response_body["error"]["message"] = "Feature '" + feature + "' is not supported by this model";
    res.set_content(response_body.dump(), "application/json");
    res.status = 400;
  }

  // Get model from name
  Model& get_model_from_name(std::vector<Model*>& models, std::string name) {
    for (auto& model : models) {
      if (model->GetName() == name) {
        return *model;
      }
    }
    throw std::runtime_error("Model not found");
  }

  // Check if model exists and return error in httplib response
  bool check_model_exists(std::vector<Model*>& models, std::string name, httplib::Response& res) {
    try {
      get_model_from_name(models, name);
    } catch (std::runtime_error& e) {
      json response_body;
      response_body["error"]["type"] = "ModelNotFound";
      response_body["error"]["message"] = "Model '" + name + "' not supported by this server!";
      res.set_content(response_body.dump(), "application/json");
      res.status = 400;
      return false;
    }
    return true;
  }

  // Provides access to a model via network
  void serveModels(std::vector<Model*> models, std::string host, int port) {

    httplib::Server svr;
    std::mutex model_mutex; // Ensure the underlying model is only called sequentially

    svr.Post("/Evaluate", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      if (!model.SupportsEvaluate()) {
        write_unsupported_feature_response(res, "Evaluate");
        return;
      }

      std::vector<std::vector<double>> inputs(request_body["input"].size());
      for (std::size_t i = 0; i < inputs.size(); i++) {
        inputs[i] = request_body["input"][i].get<std::vector<double>>();
      }

      json empty_default_config;
      json config_json = request_body.value("config", empty_default_config);

      if (!check_input_sizes(inputs, config_json, model, res))
        return;

      const std::lock_guard<std::mutex> model_lock(model_mutex);
      std::vector<std::vector<double>> outputs = model.Evaluate(inputs, config_json);

      if (!check_output_sizes(outputs, config_json, model, res))
        return;

      json response_body;
      for (std::size_t i = 0; i < outputs.size(); i++) {
        response_body["output"][i] = outputs[i];
      }

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Post("/Gradient", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      if (!model.SupportsGradient()) {
        write_unsupported_feature_response(res, "Gradient");
        return;
      }

      unsigned int inWrt = request_body.at("inWrt");
      unsigned int outWrt = request_body.at("outWrt");

      std::vector<std::vector<double>> inputs(request_body["input"].size());
      for (std::size_t i = 0; i < inputs.size(); i++) {
        inputs[i] = request_body["input"][i].get<std::vector<double>>();
      }

      std::vector<double> sens = request_body.at("sens");

      json empty_default_config;
      json config_json = request_body.value("config", empty_default_config);

      if (!check_input_wrt(inWrt, config_json, model, res))
        return;
      if (!check_output_wrt(outWrt, config_json, model, res))
        return;
      if (!check_input_sizes(inputs, config_json, model, res))
        return;
      if (!check_sensitivity_size(sens, outWrt, config_json, model, res))
        return;

      const std::lock_guard<std::mutex> model_lock(model_mutex);
      std::vector<double> gradient = model.Gradient(outWrt, inWrt, inputs, sens, config_json);

      json response_body;
      response_body["output"] = gradient;

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Post("/ApplyJacobian", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      if (!model.SupportsApplyJacobian()) {
        write_unsupported_feature_response(res, "ApplyJacobian");
        return;
      }

      unsigned int inWrt = request_body.at("inWrt");
      unsigned int outWrt = request_body.at("outWrt");

      std::vector<std::vector<double>> inputs(request_body["input"].size());
      for (std::size_t i = 0; i < inputs.size(); i++) {
        inputs[i] = request_body["input"][i].get<std::vector<double>>();
      }

      std::vector<double> vec = request_body.at("vec");

      json empty_default_config;
      json config_json = request_body.value("config", empty_default_config);

      if (!check_input_wrt(inWrt, config_json, model, res))
        return;
      if (!check_output_wrt(outWrt, config_json, model, res))
        return;
      if (!check_input_sizes(inputs, config_json, model, res))
        return;
      if (!check_vector_size(vec, inWrt, config_json, model, res))
        return;

      const std::lock_guard<std::mutex> model_lock(model_mutex);
      std::vector<double> jacobian_action = model.ApplyJacobian(outWrt, inWrt, inputs, vec, config_json);

      json response_body;
      response_body["output"] = jacobian_action;

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Post("/ApplyHessian", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      if (!model.SupportsApplyHessian()) {
        write_unsupported_feature_response(res, "ApplyHessian");
        return;
      }

      unsigned int outWrt = request_body.at("outWrt");
      unsigned int inWrt1 = request_body.at("inWrt1");
      unsigned int inWrt2 = request_body.at("inWrt2");

      std::vector<std::vector<double>> inputs(request_body["input"].size());
      for (std::size_t i = 0; i < inputs.size(); i++) {
        inputs[i] = request_body["input"][i].get<std::vector<double>>();
      }

      std::vector<double> sens = request_body.at("sens");
      std::vector<double> vec = request_body.at("vec");

      json empty_default_config;
      json config_json = request_body.value("config", empty_default_config);

      if (!check_input_wrt(inWrt1, config_json, model, res))
        return;
      if (!check_input_wrt(inWrt2, config_json, model, res))
        return;
      if (!check_output_wrt(outWrt, config_json, model, res))
        return;
      if (!check_input_sizes(inputs, config_json, model, res))
        return;
      if (!check_sensitivity_size(sens, outWrt, config_json, model, res))
        return;

      const std::lock_guard<std::mutex> model_lock(model_mutex);
      std::vector<double> hessian_action = model.ApplyHessian(outWrt, inWrt1, inWrt2, inputs, sens, vec, config_json);

      json response_body;
      response_body["output"] = hessian_action;

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Get("/Info", [&](const httplib::Request &, httplib::Response &res) {
      json response_body;
      response_body["protocolVersion"] = 1.0;
      std::vector<std::string> model_names;
      for (auto& model : models) {
        model_names.push_back(model->GetName());
      }
      response_body["models"] = model_names;

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Post("/ModelInfo", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      json response_body;
      response_body["support"] = {};
      response_body["support"]["Evaluate"] = model.SupportsEvaluate();
      response_body["support"]["Gradient"] = model.SupportsGradient();
      response_body["support"]["ApplyJacobian"] = model.SupportsApplyJacobian();
      response_body["support"]["ApplyHessian"] = model.SupportsApplyHessian();

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Post("/InputSizes", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      json empty_default_config;
      json config_json = request_body.value("config", empty_default_config);

      json response_body;
      response_body["inputSizes"] = model.GetInputSizes(config_json);

      res.set_content(response_body.dump(), "application/json");
    });

    svr.Post("/OutputSizes", [&](const httplib::Request &req, httplib::Response &res) {
      json request_body = json::parse(req.body);
      if (!check_model_exists(models, request_body["name"], res))
        return;
      Model& model = get_model_from_name(models, request_body["name"]);

      json empty_default_config;
      json config_json = request_body.value("config", empty_default_config);

      json response_body;
      response_body["outputSizes"] = model.GetOutputSizes(config_json);

      res.set_content(response_body.dump(), "application/json");
    });

    std::cout << "Listening on port " << port << "..." << std::endl;
    svr.listen(host.c_str(), port);
    std::cout << "Quit" << std::endl;
  }

}

#endif
