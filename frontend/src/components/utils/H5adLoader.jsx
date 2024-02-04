// this component with a .h5ad file as input, and return a json file for visualization

import * as hdf5 from 'jsfive'
import React, { useRef, useState, useEffect } from 'react'

const H5adLoader = ({ file }) => {
  const reader = new FileReader()
  reader.onload = (event) => {
    var f = new hdf5.File(event.target.result, file.name)
    let g = f.get('obs')
    console.log(g)
  }
  reader.readAsArrayBuffer(file)
}
