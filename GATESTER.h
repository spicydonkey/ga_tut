// A library of Mathematical functions to test a genetic algorithm
// <cmath> library must be imported

float schwefel (float *argv, int dim)
{
	float value = 418.9829f*(float)dim;
	for (int i=0; i<dim; i++)
	{
		value -= argv[i]*sin(sqrt(fabs(argv[i])));
	}
	return value;
}
