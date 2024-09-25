import React, { useState } from 'react'
import { UserOutlined } from '@ant-design/icons'
import { AutoComplete, Input, Col, Row, Popover, ConfigProvider } from 'antd'
import { useNavigate, Link } from 'react-router-dom'
import Loading from '../charts/Loading'
import axios from 'axios'
import TextCollapse from './TextCollapse'
import { useMutation } from 'react-query'
const { enterLoading, quitLoading } = Loading
const renderTitle = (title) => (
  <span>
    {title}
    <Link
      style={{
        float: 'right',
      }}
      to="/database"
    >
      more
    </Link>
  </span>
)
const renderItem = (title, key) => ({
  value: key,
  label: (
    <div
      style={{
        display: 'flex',
        justifyContent: 'space-between',
      }}
    >
      <Popover content={title} zIndex={2000}>
        {title.length > 90 ? `${title.substring(0, 90)}...` : title}
      </Popover>
      <span>
        <UserOutlined /> {key}
      </span>
    </div>
  ),
})
const options = [
  {
    label: renderTitle('Spatial Transcriptomics Studies'),
    options: [renderItem('Spatial transcriptomics map of the embryonic mouse brain – a tool to explore neurogenesis', 'STDS0000235'),
    renderItem('Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex inflammatory cell network in atopic dermatitis [ST]', 'STDS0000212'),
    renderItem("Spatial transcriptomics of de novo NEPC and ARPC", 'STDS0000227'),
    renderItem("Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer", 'STDS0000223'),
    renderItem("Spatial resolution of cellular senescence dynamics in colorectal liver metastasis", 'STDS0000219'),
    renderItem("Spatial transcriptomics of adenoid cystic carcinoma of the lacrimal gland", 'STDS0000221'),
    renderItem("Spatial transcriptomics (Visium, 10x Genomics) data of Duchenne mouse models", 'STDS0000204'),
    renderItem("Discovering Haematoma-Stimulated Circuits for Secondary Brain Injury after Intraventricular Haemorrhage by Spatial Transcriptome Analysis", 'STDS0000201'),
    renderItem("Moxibustion Improves Hypothalamus Aqp4 Polarization in APP/PS1 Mice: Evidence from Spatial Transcriptomics", 'STDS0000189'),
    renderItem("A topographic atlas defines developmental origins of cell heterogeneity in the human embryonic lung [ST]", "STDS0000175"),
    renderItem("Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment", 'STDS0000153'),
    renderItem("Intratumor heterogeneity and T cell exhaustion in primary CNS lymphoma", 'STDS0000152'),
    renderItem("Spatially resolved transcriptomics revised human distal lung epithelial hierarchy [spatial transcriptomics]", 'STDS0000114'),
    renderItem("SAW: An efficient and accurate data analysis workflow for Stereo-seq spatial transcriptomics]", 'STDS0000234')
    ],
  },
]


const Search = ({ setLocState, setCompLoad, setSelectedKey }) => {
  const [loadText, setLoadText] = useState("Browse")
  const navigate = useNavigate()
  const [loading, setLoading] = useState([])
  const [open, setOpen] = useState(false)
  const [status, setStatus] = useState("Success")
  const fetch = useMutation({
    mutationKey: ['example'],
    mutationFn: (e) => {
      axios({
        method: 'GET',
        url: '/api/example',
        params: {
          id: e
        },
      }).then((response) => {
        let dataCol = response.data.attributes
        let spitem = response.data.data[0]
        let values = Object.fromEntries(
          dataCol.map((k, i) => [k, spitem[i]])
        )
        let scitem = response.data.data[1]
        let state = {
          st: values
        }
        if (typeof scitem !== 'undefined') {
          let scvalues = Object.fromEntries(
            dataCol.map((k, i) => [k, scitem[i]])
          )
          state['sc'] = scvalues
        }
        quitLoading(0, setLoading)
        setLoadText("Browse")
        setLocState(state)
        setCompLoad({  // reset ComponentLoad
          "LassoView": false,
          "LayerView": false,
          "PairView": false,
          "MarkerTable": false,
          "NetWorkRelation": false,
          "CPDBHeatmap": false,
        })
        setSelectedKey("MetaInfo")
      })
    }
  })
  return (
    <ConfigProvider theme={{
      components: {
        Popover: {
          colorBgElevated: '#ddf7d4',
        }
      }
    }}>
      <AutoComplete
        popupClassName="certain-category-search-dropdown"
        popupMatchSelectWidth="80%"
        style={{
          width: "100%"
        }}
        options={options}
        open={open}
        onBlur={() => { setOpen(false) }}
        onSelect={() => { setOpen(false) }}
        size="large"
      >
        <Input.Search
          size="large"
          placeholder="Search for a Study to Browse."
          enterButton={loadText}
          status={status}
          loading={loading[0]}
          onClick={() => { setOpen(!open) }}
          onSearch={(e) => {
            if (e.length > 0) {
              setOpen(false)
              setStatus("Success")
              enterLoading(0, setLoading)
              setLoadText("Fetching data...")
              fetch.mutate(e)
            } else {
              setStatus("error")
            }

          }} />
      </AutoComplete>
    </ConfigProvider>
  )
}

export default Search
