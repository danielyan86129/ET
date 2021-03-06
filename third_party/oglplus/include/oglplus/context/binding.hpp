/**
 *  @file oglplus/context/binding.hpp
 *  @brief Wrappers for object binding operations
 *
 *  @author Matus Chochlik
 *
 *  Copyright 2010-2014 Matus Chochlik. Distributed under the Boost
 *  Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#pragma once
#ifndef OGLPLUS_CONTEXT_BINDING_1404142131_HPP
#define OGLPLUS_CONTEXT_BINDING_1404142131_HPP

#include <oglplus/config_compiler.hpp>
#include <oglplus/bound.hpp>

namespace oglplus {
namespace context {

/// Wrapper for object binding operations
/**
 *  @ingroup ogl_context
 */
class ObjectBinding
{
public:
	/// Binds the specified object to the specified target
	/**
	 *  Equivalent to:
	 *  @code
	 *  object.Bind(target);
	 *  @endcode
	 */
	template <typename Object>
	static void Bind(typename Object::Target target, const Object& object)
	{
		object.Bind(target);
	}

	/// Binds the specified object to the specified indexed target
	/**
	 *  Equivalent to:
	 *  @code
	 *  object.Bind(target, index);
	 *  @endcode
	 */
	template <typename Object>
	static void Bind(
		typename Object::IndexedTarget target,
		GLuint index,
		const Object& object
	)
	{
		object.Bind(target, index);
	}

	/// Binds the specified object to the appropriate binding point
	/**
	 *  Equivalent to:
	 *  @code
	 *  object.Bind();
	 *  @endcode
	 */
	template <typename Object>
	static void Bind(const Object& object)
	{
		object.Bind();
	}

	/// Uses (makes current) the specified object
	/**
	 *  Equivalent to:
	 *  @code
	 *  object.Use();
	 *  @endcode
	 */
	template <typename Object>
	static void Use(const Object& object)
	{
		object.Use();
	}

	/// Returns a managed reference to the object currently bound to target
	template <typename ObjectTarget>
	static oglplus::Bound<typename ObjectTargetOps<ObjectTarget>::Type>
	Current(ObjectTarget target)
	{
		typedef typename ObjectTargetOps<ObjectTarget>::Type ObjectOps;

		GLuint name = BindingQuery<ObjectOps>::QueryBinding(target);
		return oglplus::Bound<ObjectOps>(
			Managed<ObjectOps>(name),
			target
		);
	}

	/// Binds the object to the specified target, returns a managed reference
	template <typename Object>
	static oglplus::Bound<
		typename ObjectTargetOps<typename Object::Target>::Type
	> Bound(typename Object::Target target, const Object& object)
	{
		object.Bind(target);
		return Current(target);
	}

	/// Returns a managed reference to object currently bound to indexed target
	template <typename ObjectTarget>
	static Managed<typename ObjectTargetOps<ObjectTarget>::Type>
	Current(ObjectTarget target, GLuint index)
	{
		typedef typename ObjectTargetOps<ObjectTarget>::Type ObjectOps;

		GLuint name = BindingQuery<ObjectOps>::QueryBinding(target, index);
		return Managed<ObjectOps>(name);
	}

	/// Returns a managed reference to the currently bound Object
	template <typename Object>
	static Managed<typename ObjectBaseOps<Object>::Type> Current(void)
	{
		typedef typename ObjectBaseOps<Object>::Type ObjectOps;
		typename Object::Target target = typename Object::Target();
		GLuint name = BindingQuery<ObjectOps>::QueryBinding(target);
		return Managed<ObjectOps>(name);
	}
};

} // namespace context
} // namespace oglplus

#endif // include guard
