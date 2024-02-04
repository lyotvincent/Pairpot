import React,{Component} from 'react';
import * as echarts from 'echarts';
import "../theme/dark"
import "../theme/vintage"

class SankeyPlot extends Component {

    static defaultProps = {
        theme: 'dark',
        title: 'SankeyPlot'
    }

    constructor(props) {
        super(props);
        console.log(props);
    }

    componentDidMount () {
        setTimeout(() => {
        this.getOption()
        })
    }

    getOption = () => {
        var chartDom = document.getElementById(this.props.id);
        var myChart = echarts.init(chartDom, this.props.theme);
        var option;
        
        option = {
            title: {
                text: this.props.title,
                left: 'center'
            },
            grid: {
                height: '70%',
                top:'15%'
            },
            series: {
                type: 'sankey',
                layout: 'none',
                emphasis: {
                focus: 'adjacency'
                },
                data: [
                {
                    name: 'a'
                },
                {
                    name: 'b'
                },
                {
                    name: 'a1'
                },
                {
                    name: 'a2'
                },
                {
                    name: 'b1'
                },
                {
                    name: 'c'
                }
                ],
                links: [
                {
                    source: 'a',
                    target: 'a1',
                    value: 5
                },
                {
                    source: 'a',
                    target: 'a2',
                    value: 3
                },
                {
                    source: 'b',
                    target: 'b1',
                    value: 8
                },
                {
                    source: 'a',
                    target: 'b1',
                    value: 3
                },
                {
                    source: 'b1',
                    target: 'a1',
                    value: 1
                },
                {
                    source: 'b1',
                    target: 'c',
                    value: 2
                }
                ]
            }
        };
        
        option && myChart.setOption(option);
    }

    render () {
        return (
            <div id={this.props.id} style={{ height: '500px', width: "550px" }}>
            </div>
        )
    }
}

export default SankeyPlot;