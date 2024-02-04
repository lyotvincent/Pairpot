import React, { useRef, useEffect } from 'react'
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

const ScHeatmap = ({ theme, title, height, width, margin }) => {
  const chartRef = useRef(null) // get current DOM container

  useEffect(() => {
    var myChart = echarts.init(chartRef.current, theme) //init the echart container
    var _data = require('../../assets/data/Mouse-Brain-heatmap-sc.json')
    var _den = require('../../assets/data/Mouse-Brain-dendrogram-sc.json')

    var den_data = []
    for (let i = 0; i < _den.icoord.length; i++) {
      den_data.push({
        coords: _den.icoord[i].map((item, index) => {
          return [item, _den.dcoord[i][index]]
        }),
      })
    }

    var cells = _den.ivl
    var genes_unsorted = _data.map((item) => {
      return item.gene
    })
    var genes = []
    var top = 5
    for (var i = 0; i < cells.length; i++) {
      genes = [
        ...genes,
        ...genes_unsorted.slice(
          _den.leaves[i] * top,
          (_den.leaves[i] + 1) * top
        ),
      ]
    }
    genes = genes.reverse()

    var expr = _data.map((item) => {
      return item.expr
    })

    myChart.setOption({
      tooltip: {},
      title: [
        {
          text: title,
          left: '40%',
          textStyle: {
            fontSize: 24,
          },
        },
      ],
      grid: [
        {
          top: '12%',
          left: '15%',
          right: '21%',
          bottom: '13%',
          axisLine: {
            lineStyle: {
              color: '#fff',
            },
          },
          axisPointer: {
            lineStyle: {
              color: '#ffbd67',
            },
          },
        },
        {
          top: '4%',
          left: '15%',
          right: '21%',
          height: '8%',
          show: false,
        },
      ],
      yAxis: [
        {
          gridIndex: 0,
          name: 'Top-5 Markers',
          nameLocation: 'middle',
          nameGap: 70,
          nameRotate: 90,
          nameTextStyle: {
            fontWeight: 'bold',
            fontSize: 16,
            padding: 1,
          },
          type: 'category',
          data: genes,
          splitArea: {
            show: true,
            interval: 4,
          },
          axisTick: {
            alignWithLabel: true,
          },
          axisLabel: {
            interval: 0,
            rotate: 0,
            margin: 10,
          },
        },
        {
          gridIndex: 1,
          type: 'value',
          min: 0,
          max: 2,
          show: false,
          axisTick: {
            show: false,
          },
          axisLabel: {
            show: false,
          },
          axisLine: {
            show: false,
          },
        },
      ],
      xAxis: [
        {
          gridIndex: 0,
          name: 'Cell-types',
          nameLocation: 'middle',
          nameGap: 70,
          nameRotate: 0,
          nameTextStyle: {
            fontWeight: 'bold',
            fontSize: 16,
            padding: 1,
          },
          type: 'category',
          data: cells,
          splitArea: {
            show: true,
          },
          axisTick: {
            alignWithLabel: true,
          },
          axisLabel: {
            interval: 0,
            fontSize: 14,
            rotate: -20,
            margin: 12,
          },
        },
        {
          gridIndex: 1,
          type: 'value',
          min: 0,
          max: 70,
          show: false,
          axisTick: {
            show: false,
          },
          axisLabel: {
            show: false,
          },
          axisLine: {
            show: false,
          },
        },
      ],
      dataset: {
        dimensions: ['cell', 'gene', 'expr', 'frac'],
        source: _data.map((item) => {
          return [item.cell, item.gene, item.expr, item.frac]
        }),
        xAxis: 0,
        yAxis: 1,
      },
      series: [
        {
          type: 'scatter',
          xAxisIndex: 0,
          yAxisIndex: 0,
          encode: {
            x: 'cell',
            y: 'gene',
            tooltip: [0, 1, 2, 3],
          },
          symbolSize: (val) => {
            return val[3] * 18
          },
          categories: cells,
          emphasis: {
            itemStyle: {
              borderColor: '#fff',
              borderWidth: 1,
            },
          },
        },
        {
          type: 'lines',
          xAxisIndex: 1,
          yAxisIndex: 1,
          coordinateSystem: 'cartesian2d',
          polyline: true,
          data: den_data,
          lineStyle: { color: '#fff', width: 2 },
          drawPoint: true,

          silent: false,
        },
      ],
      visualMap: [
        {
          left: 'right',
          top: '23%',
          orient: 'vertical',
          calculable: true,
          dimension: 'frac',
          itemWidth: 20,
          itemHeight: 120,
          inRange: {
            // 设置映射的颜色范围为蓝色到红色之间
            symbolSize: [3, 18],
            color: ['#e3e3e3', '#111111'],
          },
          text: ['Fraction of cells\n in group(%)'],
          textGap: 20,
          textStyle: {
            fontSize: 14,
          },
          padding: 8,
          precision: 1,
          label: {
            formatter: '{value} %',
          },
          min: 0,
          max: 1,
          seriesIndex: 0,
          handleIcon: 'circle',
          handleSize: '100%',
        },
        {
          left: 'right',
          bottom: '23%',
          orient: 'vertical',
          calculable: true,
          precision: 1,
          min: 0,
          max: Math.max.apply(null, expr),
          itemHeight: 120,
          dimension: 'expr',
          seriesIndex: 0,
          text: ['Mean Expression\nin group'],
          textGap: 13,
          textStyle: {
            fontSize: 14,
          },
          indicatorIcon: 'roundRect',
          indicatorSize: '60%',
        },
      ],
    })
  }, [theme, title])

  return (
    <div>
      <div
        ref={chartRef}
        className="chart"
        //the target DOM container needs height and width
        style={{ height: height, width: width, margin: margin }}></div>
    </div>
  )
}

ScHeatmap.defaultProps = {
  theme: 'dark',
  title: 'ScHeatmap',
  height: '50rem',
  width: '45rem',
  margin: '2rem',
}

ScHeatmap.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default ScHeatmap
