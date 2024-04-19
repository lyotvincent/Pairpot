// Correlations plots among cell-types via CellPhoneDB, including Heatmap, SankeyPlot, GraphPlot
import React, { useRef, useEffect, useState } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
  TooltipComponent,
  VisualMapComponent,
} from 'echarts/components'
import { LineChart, LinesChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import { ContactsOutlined } from '@ant-design/icons'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  LineChart,
  CanvasRenderer,
  UniversalTransition,
  VisualMapComponent,
  LinesChart,
])

const CpdbCorr = ({ theme, height, width, margin }) => {

  const [isInit, setInit] = useState(false) // whether echart object is inited
  const chartRef = useRef(null) // current DOM container
  const [_data, _setData] = useState(require('../../assets/data/significant_means.json'))

  const CpdbTsvLoader = (file) => {
    const reader = new FileReader()
    reader.onload = (e) => {
      let fileContent = e.target.result
      let rows = fileContent.split('\n').map(row => row.split('\t'))
      _setData(rows)
    }
    reader.readAsArrayBuffer(file)
  }

  useEffect(() => {
    var chart = echarts.init(chartRef.current, theme) //init the echart container
    let _dims = [...Object.keys(_data)]
    let _typeDims = _dims.slice(14)
    let _type = new Set(_typeDims.map((item) => item.split('|')[0]))
    _type = Array.from(_type)
    console.log(_type)
    let _sum = _typeDims.map((item) => _data[item].filter((i) => i !== null).length)
    let _sum2D = _type.map((_, index) => _sum.slice(index * _type.length, (index + 1) * _type.length))
    let _typeCum = _type.map((_, index) => _sum2D[index].reduce((acc, cur) => acc + cur, 0))
    let _typeCumCol = _type.map((_, index) => _sum2D.map((item) => item[index]).reduce((acc, cur) => acc + cur, 0))
    _typeCum = _typeCum.map((item, index) => _typeCumCol[index] + item)
    let _typeIdx = _typeCum.map((_, index) => index)
    _typeIdx = _typeIdx.sort((a, b) => _typeCum[a] - _typeCum[b])
    let elem = _sum.map((_, index) => {
      let x = Math.floor(index / _type.length)
      let y = index % _type.length
      let val = _sum2D[_typeIdx[x]][_typeIdx[y]]
      return ([x, y, val])
    })
    // prettier-ignore
    const hours = _typeIdx.map((item) => _type[item])
    // prettier-ignore
    const days = _typeIdx.map((item) => _type[item])
    // prettier-ignore
    const data = elem
    chart.setOption({
      title: {
        text: "Heatmap",
        left: 'center',
      },
      tooltip: {
        position: 'top',
      },
      grid: {
        top: '10%',
        left: '10%',
        right: '30%'
      },
      xAxis: {
        type: 'category',
        axisTick: {
          interval: 1
        },
        data: hours,
        splitArea: {
          show: true,
        },
        axisLabel: {
          rotate: 90
        }
      },
      yAxis: {
        type: 'category',
        position: 'right',
        inverse: true,
        data: days,
        splitArea: {
          show: true,
        },
        axisLabel: {
          rotate: 0
        }
      },
      visualMap: {
        min: 0,
        max: 20,
        calculable: true,
        orient: 'vertical',
        left: '0%',
        bottom: 'center',
        inRange: {
          color: ["#104e8b", "#ffdab9", "#8b0a50"],
        },
      },
      series: [
        {
          name: 'Punch Card',
          type: 'heatmap',
          data: data,
          label: {
            show: true,
          },
          emphasis: {
            itemStyle: {
              shadowBlur: 10,
              shadowColor: 'rgba(0, 0, 0, 0.5)',
            },
          },
        },
      ],
    })
  }, [theme])


  return (
    <div
      ref={chartRef}
      className="chart"
      //the target DOM container needs height and width
      style={{ height: height, width: width, margin: margin }}></div>
  )
}

CpdbCorr.defaultProps = {
  theme: 'dark',
  height: '31rem',
  width: '40rem',
  margin: '1rem',
}

CpdbCorr.propTypes = {
  theme: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default CpdbCorr
