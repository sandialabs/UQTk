import umbridge
import math

class TestModel(umbridge.Model):

    def __init__(self):
        super().__init__("forward")

    def get_input_sizes(self, config):
        return [3]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        x = parameters[0]

        lik = math.exp(-.5*((x[0]-0.5)*(x[0]-0.5)/.01 + (x[1]-0.5)*(x[1]-0.5)/.01 + (x[2]-0.5)*(x[2]-0.5)/.01)) \
            + math.exp(-.5*((x[0]+0.5)*(x[0]+0.5)/.01 + (x[1]+0.5)*(x[1]+0.5)/.01 + (x[2]+0.5)*(x[2]+0.5)/.01)) \
            + math.exp(-.5*((x[0]-0.2)*(x[0]-0.2)/.01 + (x[1]+0.2)*(x[1]+0.2)/.01 + (x[2]+0.2)*(x[2]+0.2)/.01))

        if (lik == 0):
            return [[-1e100]]
        return [[math.log(lik)]]

    def supports_evaluate(self):
        return True

testmodel = TestModel()

umbridge.serve_models([testmodel], 4242)
