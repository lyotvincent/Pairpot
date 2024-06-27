const setEmptyAxis = (gridIndex) => ({
  xAxis: [{
    id: gridIndex,
    gridIndex: gridIndex,
    name: gridIndex,
    nameLocation: 'middle',
    nameGap: 23,
    nameTextStyle: {
      fontSize: 16,
    },
    axisLine: {
      onZero: false,
    },
    axisLabel: {
      fontSize: 14,
    },
    position: 'bottom',
  }],
  yAxis:[{
    id: gridIndex,
    gridIndex: gridIndex,
    name: gridIndex,
    nameTextStyle: {
      fontSize: 16,
    },
    axisLine: {
      onZero: false,
    },
    axisLabel: {
      fontSize: 14,
    },
    position: 'left',
  }],
})

const setEmptyAxis3D = (gridIndex) =>({
  xAxis3D: [{
    grid3DIndex: gridIndex,
    name: gridIndex,
    nameLocation: 'middle',
  }
  ],
  yAxis3D: [{
    grid3DIndex: gridIndex,
    name: gridIndex,
    nameLocation: 'middle',
  }],
  zAxis3D: [{
    gridIndex: gridIndex,
    name: gridIndex,
    nameLocation: 'top',
  }],
})

const Axis = {
  setEmptyAxis,
  setEmptyAxis3D,
}

export default Axis