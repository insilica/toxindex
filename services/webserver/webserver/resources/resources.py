import os
import dotenv
import pathlib


package_directory = os.path.dirname(os.path.abspath(__file__))

def path(path):
  return os.path.join(package_directory,path)

def get_env(key=None):
  env = dotenv.dotenv_values(path('.env'))
  return env if key is None else env[key]

# TODO should replace this with flask.url_for?
def abs_route(path):
  return f"{get_env('ROOT_URL')}/{path}"