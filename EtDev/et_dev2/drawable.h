// base class for renderable
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DRAWABLE_H
#define DRAWABLE_H

#include "all_gl.h"

#include <cassert>
#include <QMouseEvent>

/// Base class for drawable class that can be attached to the QGLWidget
class Drawable
{
public:
	//Drawable();
	virtual ~Drawable(void) {}

	/*/// prepare for future rendering the given mesh object
	virtual void setMeshToDraw(std::shared_ptr<MyMesh> _m) = 0
	{

	}*/

	virtual void reset() = 0
	{
		assert(!"Drawable::reset() cannot be called directly!");
	}

	virtual void setLightingEnabled(bool _is_enabled) 
	{

	}

	virtual bool getLightingEnabled()
	{
		return true;
	}

	virtual void setPerVertColor(float* _color_data, int _n) = 0
	{
		assert(!"Drawable::setPerVertColor() cannot be called directly!");
	}

	virtual void setTransparencyEnabled(const bool _enabled)
	{
		
	}

	virtual bool getTransparencyEnabled() 
	{
		return true;
	}

	/// Rendering procedure with simple timing
	virtual void render(double /*time*/) = 0
	{
		assert(!"Drawable::render() cannot be called directly!");
	}

	/// Reshape event handler
	virtual void reshape(GLuint _w, GLuint _h) = 0
	{
		assert(!"Drawable::Reshape cannot be called directly!");
	}

	/// clear the bound buffer data
	virtual void clear()
	{

	}

	/// respond to mouse movement
	virtual void MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h)
	{
		// TODO
	}

	virtual void MousePress(QMouseEvent* _evnt)
	{
		// TODO
	}

	virtual void MouseRelease(QMouseEvent* _evnt)
	{
		// TODO
	}
};


#endif