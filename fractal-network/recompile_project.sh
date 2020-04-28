BUILD_DIR="build"
BUILD_TYPE="Release"

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "Directory \"$BUILD_DIR\" does not exist. Creating."
  mkdir "$BUILD_DIR"
fi

cd "${BUILD_DIR}"; cmake -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" ..; make
