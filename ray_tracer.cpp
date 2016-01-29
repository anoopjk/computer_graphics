/**
 * ray_tracer.cpp
 * CS230 Assignment 2, Winter 2012
 * -------------------------------
 * Implement ray tracer here.
 */
#include <algorithm> // include algorithm for max
#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))


#include "ray_tracer.h"

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}
///////////////////////////////////////////////////////////////////////////
// Lambertian model
//Ld = kd*I*max(0, n · l)
Vector_3D<double> Lambertian_shading(const Ray& ray,
									 const Light* light,
									 const Vector_3D<double>& color_diffuse,
									 const Vector_3D<double>& intersection_point,
									 const Vector_3D<double>& same_side_normal)
{
	
		Vector_3D<double> light_position = light->position;
		Vector_3D<double> light_ray = light_position - intersection_point;
		light_ray.Normalize();
		double ndotl = Vector_3D<double>::Dot_Product(light_ray, same_side_normal); // n.l
		double maximum = std::max(0.0,ndotl);
    	Vector_3D<double> Intensity_em = light->Emitted_Light(ray);
		//Ld = kd*I*max(0, n · l)
		Vector_3D<double> color;
		color.x += color_diffuse.x*Intensity_em.x*maximum; //kd*I*max(0, n · l)
		color.y += color_diffuse.y*Intensity_em.y*maximum;
		color.z += color_diffuse.z*Intensity_em.z*maximum;
	
		return color;
	}
	
// specular model
//Ls = ks*I*max(0, n · h)^p ,
Vector_3D<double> Specular_shading(const Ray& ray,
								   const Light* light,
								   const Vector_3D<double>& color_specular,
								   const Vector_3D<double>& intersection_point,
								   const Vector_3D<double>& same_side_normal, 
								   const double specular_power)
{
	
	Vector_3D<double> light_ray = light->position - intersection_point;
	Vector_3D<double> view_ray = ray.endpoint - intersection_point;
	view_ray.Normalize();
	Vector_3D<double> h = (view_ray + light_ray);
	h.Normalize();
    
	double ndoth = Vector_3D<double>::Dot_Product(h, same_side_normal);
	double maximum = std::max(0.0, ndoth);
	//Ls = ks*I*max(0, n · h)^p ,
	Vector_3D<double> color;
	Vector_3D<double> Intensity_em = light->Emitted_Light(ray);
	color.x += color_specular.x*Intensity_em.x*(pow(maximum,specular_power)); //ks*I*max(0, n · h)^p
	color.y += color_specular.y*Intensity_em.y*(pow(maximum, specular_power));
	color.z += color_specular.z*Intensity_em.z*(pow(maximum, specular_power));

	return color;
	
	
}
//------------------------------------------------------------------------------------
// Shader
//----------------------------------------------------------------------------------
// sum of ambient, diffuse and specular color models
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,
              const Object& intersection_object,
			  const Vector_3D<double>& intersection_point,
			  const Vector_3D<double>& same_side_normal) const
{
    
	Vector_3D<double> color;
	//Vector_3D<double>lambert_shading;
	//Vector_3D<double> specular_shading;
	Vector_3D<double> ambient_shading;
	Light* light;
// for all the lights in world find lambertian and specular and sum them (superposition) at last add ambient_color
	for (unsigned int i = 0; i < this->world.lights.size(); ++i)
	{
	    light = this->world.lights.at(i);
		Vector_3D<double> light_position = light->position; // endpoint of light
		Vector_3D<double> light_direction = (light_position - intersection_point); // direction of light ray
		light_direction.Normalize();
		Vector_3D<double> p = intersection_point + light_direction*Object::small_t; // numerical imprecision causing the ray to hit the surface p is on.
		Ray light_ray(p, light_direction); // generate light rays from the points on the ground plane
		bool shadow_ray = false;
		for (unsigned int j = 0; j < world.objects.size(); ++j)
		{
			if (world.objects.at(j)->Intersection(light_ray)) // the object has intersection point with the generated ray p then it is a shadowray
				shadow_ray = true;
		}
		//cout << intersection_object.Intersection(light_ray) << endl;
		//cout << intersection_object.Normal(intersection_point) << endl;
		if (!shadow_ray || !(this->world.enable_shadows)){
			Vector_3D<double> lambert_shading = Lambertian_shading(ray, light, this->color_diffuse, intersection_point, same_side_normal);
			Vector_3D<double> specular_shading = Specular_shading(ray, light, this->color_specular, intersection_point, same_side_normal, this->specular_power);
			color += lambert_shading + specular_shading;
		}
		Vector_3D<double>Intensity_a = light->Emitted_Light(ray);
		ambient_shading += Intensity_a*(this->color_ambient);
	}
	// ambient color is added irrespective of shadows
	color += ambient_shading*0.25 ;    // ambient color is uniform. reducing it to 25% .
    	
    return color;
}
//---------------------------------------------------------------------------------------------
//Reflective shader
Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,
              const Object& intersection_object,
			  const Vector_3D<double>& intersection_point,
			  const Vector_3D<double>& same_side_normal) const
{
	Vector_3D<double> color;
	Vector_3D<double> reflect_color;
	color = Phong_Shader::Shade_Surface(ray, intersection_object, intersection_point, same_side_normal);
	if (reflectivity != 0)
	{
		 Ray reflect_ray;
		//while (reflect_ray.recursion_depth < this->world.recursion_depth_limit) {
			reflect_ray.endpoint = intersection_point;
			//r = d − 2(d · n)n  
			reflect_ray.direction = ray.direction - same_side_normal*(Vector_3D<double>::Dot_Product(ray.direction, same_side_normal)) * 2;
			reflect_ray.direction.Normalize();
			reflect_ray.t_max = FLT_MAX;  //s ∈ [ small_t, ∞)  if we set this value to small_t, there won't be any reflections
			reflect_color += world.Cast_Ray(reflect_ray, ray);
			reflect_ray.recursion_depth += 1;
		//}
	}

		//cout << reflect_ray.recursion_depth << endl;
    //color c = c + km*raycolor(p+ sr, e	,∞)
	color = color + reflect_color*reflectivity;
	return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,
              const Object& intersection_object,
			  const Vector_3D<double>& intersection_point,
			  const Vector_3D<double>& same_side_normal) const
{
    return color;
}
//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::
Intersection(Ray& ray) const
{
	//sphere -ray intersection is given by
	// (d · d)t^2 + 2*d·(e − c)*t + (e − c)·(e − c) − R^2 = 0.
	double A = 1;
	double B = 2 * (Vector_3D<double>::Dot_Product(ray.direction, (ray.endpoint - this->center)));
	double C = Vector_3D<double>::Dot_Product((ray.endpoint - this->center), (ray.endpoint - this->center)) - sqr(this->radius);
	
	
	double discriminant = sqr(B) - 4 * A*C;

	if (discriminant == 0)
	{
		double t1 = -B / (2 * A);
		if (t1 > small_t)
		{
			ray.t_max = t1;
			ray.semi_infinite = false;
			ray.current_object = this;
			return true;
		}
		return false;
	}
	else if (discriminant > 0)
	{
		double t1 = (-B + sqrt(discriminant)) / (2 * A);
		double t2 = (-B - sqrt(discriminant)) / (2 * A);
		if ((t1 < t2) && (t1> small_t))
		{
			ray.t_max = t1;
			ray.semi_infinite = false;
			ray.current_object = this;
			return true;
		}
		else if ((t1 > t2) && (t2 > small_t))
		{
			ray.t_max = t2;
			ray.semi_infinite = false;
			ray.current_object = this;
			return true;
		}
		return false;
	}
	else 
	{
		return false;
	}


}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;
    normal =  (location - center)*(1/radius); //n = (p-c)/r ;
    return normal;
}


// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{
	//t = (x1- p).N / d.N 
	//
	Vector_3D<double> ray_plane = x1 - ray.endpoint ;
	double  Nr = Vector_3D<double>::Dot_Product(ray_plane, this->normal);
	double Dr = Vector_3D<double>::Dot_Product(ray.direction, this->normal);
	//cout << Nr << endl;
    if(-Nr > small_t )
	{
		//cout << "entered" << endl;
		double t1 = Nr/Dr;
		
		if(t1 > small_t){
			ray.t_max = t1; 
			ray.semi_infinite = false; 	 
			ray.current_object = this;   
			return true;
		}	
	} 
    return false;
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
	// TODO 
	double WIDTH = this->film.pixel_grid.m;	
	double HEIGHT = this->film.pixel_grid.n;
	double left = -WIDTH/2;
	double right =-left;
	double bottom = -HEIGHT/2;
	double top = -bottom;
	

	double u = left   + ((right-left)*(pixel_index.x  + 0.5))/WIDTH;
	double v = bottom + ((top-bottom)*(pixel_index.y  + 0.5))/HEIGHT;
	

	Vector_3D<double> result = this->focal_point + 
		                       this->horizontal_vector*u*this->film.pixel_grid.dx +
		                       this->vertical_vector*v*this->film.pixel_grid.dy;                                      
									                                 
	return  result;
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{


    // Iterate over all the scene object to get intersection
   const Object* closest_obj = 0;
    double t_max = ray.t_max;
    for(unsigned int i = 0; i < this->objects.size(); ++i)
	{
		if(this->objects.at(i)->Intersection(ray)){						
			closest_obj = ray.current_object;
		} 			
	}
    return closest_obj;
}

// set up the initial view ray and call 
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    Ray ray; // TODO: set up the initial view ray here

	            ray.endpoint = this->camera.position;// perspective projection- origin is camera for all rays
				Vector_3D<double> color;
	
				ray.direction = camera.World_Position(pixel_index) - this->camera.position;
				ray.direction.Normalize();
				ray.t_max = FLT_MAX;		
				Ray dummy_root = ray;	
				color =  Cast_Ray(ray,dummy_root);
	
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}
	

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
	Vector_3D<double> color;
	const Object *closest_obj = Closest_Intersection(ray);
	   if(closest_obj != 0)
	   {
		
		Vector_3D<double> intersection_point = ray.endpoint + ray.direction*ray.t_max;
		Vector_3D<double> same_side_normal = closest_obj->Normal(intersection_point);
		
		same_side_normal.Normalize();
		color = closest_obj->material_shader->Shade_Surface(parent_ray,
			                                                *closest_obj,intersection_point,
															same_side_normal);				
	  }
	   else
	   {
		   color = background_shader->Shade_Surface(parent_ray, 
			                                         *closest_obj, 
													 Vector_3D<double>(), 
													 Vector_3D<double>());

	   }
	
    return color;
}


