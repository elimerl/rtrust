struct Uniforms {
    screen_size: vec2<f32>
};
@group(0) @binding(0) // 1.
var<uniform> uniforms: Uniforms;

var<private> a: array<vec4<f32>,3> = array(
  		vec4(-1., -1., 0., 1.),
  		vec4(3., -1., 0., 1.),
  		vec4(-1., 3., 0., 1.),
);

@vertex
fn vs_main(@builtin(vertex_index) in_vertex_index: u32) -> @builtin(position) vec4<f32> {
    return a[in_vertex_index];
}


@fragment
fn fs_main(@builtin(position) in: vec4<f32>) -> @location(0) vec4<f32> {
	let scene: Scene = Scene(array(Sphere(vec3<f32>(0.), 1.0)), 1);
	let focal_length = 1.0;
    let viewport_height = 2.0;
    let viewport_width = viewport_height * (uniforms.screen_size.x/uniforms.screen_size.y);
    let camera_center = vec3<f32>(0., 0., 0.);

    // Calculate the vectors across the horizontal and down the vertical viewport edges.
    let viewport_u = vec3<f32>(viewport_width, 0., 0.);
    let viewport_v = vec3<f32>(0., -viewport_height, 0.);

    // Calculate the horizontal and vertical delta vectors from pixel to pixel.
    let pixel_delta_u = viewport_u / uniforms.screen_size.x;
    let pixel_delta_v = viewport_v / uniforms.screen_size.y;

    // Calculate the location of the upper left pixel.
    let viewport_upper_left = camera_center
                             - vec3<f32>(0., 0., focal_length) - viewport_u/2. - viewport_v/2.;
    let pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

    let pixel_center = pixel00_loc + (in.x * pixel_delta_u) + (in.y * pixel_delta_v);
    let ray_direction = pixel_center - camera_center;
    let r = Ray(camera_center, ray_direction);
    let pixel_color = ray_color(r);
    return pow(vec4<f32>(pixel_color, 1.0), vec4(2.2));
}

struct Ray {
	origin: vec3<f32>,
	direction: vec3<f32>	
};

fn ray_at(ray: Ray, t: f32) -> vec3<f32> {
	return ray.origin + t * ray.direction;
}

fn ray_color(ray: Ray, scene: Scene) -> vec3<f32> {
	let rec = hit_scene(scene, 0.0, INFINITY);
    if (rec_is_hit(rec.t))
    {
        return 0.5 * (rec.normal + vec3<f32>(1.0));
	}

	let unit_direction = normalize(ray.direction);
	let a = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - a) * vec3<f32>(1.0, 1.0, 1.0) + a * vec3<f32>(0.5, 0.7, 1.0);
}

fn hit_sphere(sphere: Sphere, r: Ray, ray_tmin: f32, ray_tmax: f32) -> HitRec {
    let oc = r.origin - sphere.center;
    let a = length_squared(r.direction);
    let half_b = dot(oc, r.direction);
    let c = length_squared(oc) - sphere.radius*sphere.radius;
    let discriminant = half_b*half_b - a*c;

    if (discriminant < 0.0) {
        return rec_no_hit();
    } else {
		let sqrtd = sqrt(discriminant);
		var root = (-half_b - sqrtd) / a;
        if (root <= ray_tmin || ray_tmax <= root) {
            root = (-half_b + sqrtd) / a;
            if (root <= ray_tmin || ray_tmax <= root)
                {return rec_no_hit();}
        }

        let t = root;
        let p = ray_at(r, t);
        let outward_normal = (p - sphere.center) / sphere.radius;

        let rec = HitRec(p, outward_normal, t, false);

		return face_normal(rec, r, outward_normal);
    }
}


fn length_squared(v: vec3<f32>) -> f32 {
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

struct Sphere {
	center: vec3<f32>,
	radius: f32
}

struct HitRec {
	p: vec3<f32>,
	normal: vec3<f32>,
	t: f32,
	front_face: bool
}

fn face_normal(rec: HitRec, r: Ray, outward_normal: vec3<f32>) -> HitRec {
	var n = rec;
	n.front_face = dot(r.direction, outward_normal) < 0.;
	if (n.front_face) {
		n.normal = outward_normal;
	} else {
		n.normal = -outward_normal;
	}
	return n;
}

fn rec_no_hit() -> HitRec { return HitRec(vec3<f32>(0.),vec3<f32>(0.),-1.0, true); }
fn rec_is_hit(rec: HitRec) -> bool { return rec.t > 0.0; }

struct Scene {
	spheres: array<Sphere, 128>,
	spheres_count: u32
}

fn hit_scene(scene: Scene, r: Ray, ray_tmin: f32, ray_tmax: f32) -> HitRec {
	var hit_anything = false;
	var closest_so_far = ray_tmax;
	var current_rec = rec_no_hit();

	for (var i: u32 = 0u; i < scene.spheres_count; i++) {
		let sphere = scene.spheres[i];
		let rec = hit_sphere(sphere, r, ray_tmin, closest_so_far);
		if (rec_is_hit(rec)) {
			hit_anything = true;
			closest_so_far = rec.t;
		}
	}

	return current_rec;
}
const INFINITY: f32 = 1e10;
const PI: f32 = 3.1415926535897932385;

fn deg2rad(deg: f32) {
	return deg * PI / 180.0;
}