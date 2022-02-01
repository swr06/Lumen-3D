#include "TAAJitter.h"

static glm::vec2 HaltonSequenceData[64];

// inside taa 
float HaltonSequence(int Prime, int index) 
{
	float r = 0.0f;
	float f = 1.0f;
	int i = index;
	while (i > 0)
	{
		f /= Prime;
		r += f * (i % Prime);
		i = (int)std::floor(i / float(Prime));
	}
	return r;
}

void Lumen::GenerateJitterStuff()
{
	for (int i = 0; i < 64; i++) 
	{
		// 2, 3 as the primes :D
		HaltonSequenceData[i].x = HaltonSequence(2, i + 1);
		HaltonSequenceData[i].y = HaltonSequence(3, i + 1);
	}
}

glm::vec2 Lumen::GetTAAJitter(int CurrentFrame)
{
	glm::vec2 Jitter = HaltonSequenceData[CurrentFrame % 64];
	return Jitter;
}

glm::vec2 Lumen::GetTAAJitterSecondary(int CurrentFrame)
{
	glm::vec2 Jitter = HaltonSequenceData[CurrentFrame % 32];
	return Jitter;
}
