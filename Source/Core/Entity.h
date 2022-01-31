#pragma once

#include "Object.h"

#include <iostream>

namespace Lumen
{

	class Entity
	{
	public : 
		Entity(Object* object) : m_Object(object)
		{
			m_Model = glm::mat4(1.0f);
		}

		Object* const m_Object;
		glm::mat4 m_Model;
	};
}