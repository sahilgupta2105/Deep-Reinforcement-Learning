from distutils.core import setup, Extension


example_module = Extension('_physics_sim',
                           sources=['simulator_interface_wrap.cxx', 'scalar_field.cpp','circle_geometry.cpp','rigid_body.cpp','MAC.cpp','simulator.cpp'],
                           extra_compile_args=['-std=c++11', '-O3'],
							include_dirs=['/Users/sahil/local/include'],
                           library_dirs=['/Users/sahil/local/lib'],
                           libraries=['partio'])

setup (name = 'simulator',
       version = '1.0',
       author      = "gsahil@seas.upenn.edu",
       description = """Physics simulator 2d""",
       ext_modules = [sim_module],
       py_modules = ["sim"],
       )