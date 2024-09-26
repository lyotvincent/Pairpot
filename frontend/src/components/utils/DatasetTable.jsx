import {
  EyeOutlined,
  CloudDownloadOutlined,
  LinkOutlined,
  SearchOutlined,
} from '@ant-design/icons'
import React, { useEffect, useRef, useState } from 'react'
import Highlighter from 'react-highlight-words'
import { Button, Input, Space, Table, Pagination, Tooltip } from 'antd'
import { useNavigate } from 'react-router-dom'
import DatasetDescription from './DatasetDescription'
import DatasetPaired from './DatasetPaired'
import { useQuery } from 'react-query'
import axios from 'axios'
const IconTip = ({ icon, attr, onClick, placement = 'bottom' }) => (
  <Tooltip placement={placement} title={attr} arrow={false} align={'center'}>
    {React.createElement(icon, (onClick = { onClick }))}
  </Tooltip>
)

const DatasetTable = ({ src, col, oriSrc }) => {
  const navigate = useNavigate()
  const [data, setData] = useState([])
  const [descOpen, setDescOpen] = useState(false)
  const [descInfo, setDescInfo] = useState({})
  const [dataCol, setDataCol] = useState([])
  const [searchText, setSearchText] = useState('')
  const [searchedColumn, setSearchedColumn] = useState('')
  const searchInput = useRef(null)
  const handleSearch = (selectedKeys, confirm, dataIndex) => {
    confirm()
    setSearchText(selectedKeys[0])
    setSearchedColumn(dataIndex)
  }
  const handleReset = (clearFilters) => {
    clearFilters()
    setSearchText('')
  }
  const getColumnSearchProps = (dataIndex) => ({
    filterDropdown: ({
      setSelectedKeys,
      selectedKeys,
      confirm,
      clearFilters,
      close,
    }) => (
      <div
        style={{
          padding: 8,
        }}
        onKeyDown={(e) => e.stopPropagation()}>
        <Input
          ref={searchInput}
          placeholder={`Search ${dataIndex}`}
          value={selectedKeys[0]}
          onChange={(e) =>
            setSelectedKeys(e.target.value ? [e.target.value] : [])
          }
          onPressEnter={() => handleSearch(selectedKeys, confirm, dataIndex)}
          style={{
            marginBottom: 8,
            display: 'block',
          }}
        />
        <Space>
          <Button
            type="primary"
            onClick={() => handleSearch(selectedKeys, confirm, dataIndex)}
            icon={<SearchOutlined />}
            size="small"
            style={{
              width: 90,
            }}>
            Search
          </Button>
          <Button
            onClick={() => clearFilters && handleReset(clearFilters)}
            size="small"
            style={{
              width: 90,
            }}>
            Reset
          </Button>
          <Button
            type="link"
            size="small"
            onClick={() => {
              confirm({
                closeDropdown: false,
              })
              setSearchText(selectedKeys[0])
              setSearchedColumn(dataIndex)
            }}>
            Filter
          </Button>
          <Button
            type="link"
            size="small"
            onClick={() => {
              close()
            }}>
            close
          </Button>
        </Space>
      </div>
    ),
    filterIcon: (filtered) => (
      <SearchOutlined
        style={{
          color: filtered ? '#1677ff' : undefined,
        }}
      />
    ),
    onFilter: (value, record) =>
      record[dataIndex].toString().toLowerCase().includes(value.toLowerCase()),
    onFilterDropdownOpenChange: (visible) => {
      if (visible) {
        setTimeout(() => searchInput.current?.select(), 100)
      }
    },
    render: (text) =>
      searchedColumn === dataIndex ? (
        <Highlighter
          highlightStyle={{
            backgroundColor: '#ffc069',
            padding: 0,
          }}
          searchWords={[searchText]}
          autoEscape
          textToHighlight={text ? text.toString() : ''}
        />
      ) : (
        text
      ),
  })

  const initColumns = [
    Table.SELECTION_COLUMN,
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: '10%',
      sorter: (a, b) => a.key - b.key,
      sortDirections: ['descend', 'ascend'],
      fixed: 'left',
    },
    {
      title: 'Title',
      dataIndex: 'title',
      key: 'title',
      width: '40%',
      ...getColumnSearchProps('title'),
    },
    {
      title: 'Species',
      dataIndex: 'species',
      key: 'species',
      width: '20%',
    },
    {
      title: 'Tissues',
      dataIndex: 'tissues',
      key: 'tissues',
      width: '20%',
    },
    {
      title: 'Technology',
      dataIndex: 'technologies',
      key: 'technologies',
      width: '20%',
    },
    {
      title: 'PMID',
      dataIndex: 'pmid',
      key: 'pmid',
      ...getColumnSearchProps('pmid'),
      sorter: (a, b) => a.pmid.length - b.pmid.length,
      sortDirections: ['descend', 'ascend'],
      width: '20%',
    },
    {
      title: 'Action',
      key: 'operation',
      fixed: 'right',
      width: 120,
      render: () => <a>Visualization</a>,
    },
  ]

  const [dyColumns, setDyColumns] = useState(initColumns)

  const OriginResponse = useQuery({
    queryKey: ['db'],
    queryFn: () => axios.get('/api/datasets').then((response) => {
      return response.data
    }).catch((error) => {
      console.log(error)
    }),
    staleTime: Infinity,
    retry: false,
    refetchOnWindowFocus: false,
  })

  const GetscItem = (scid) => {
    return OriginResponse.data?.data.find((item) => item[1] === scid)
  }

  useEffect(() => {
    setDataCol(col)
    var attrStack = {}
    for (let i = 0; i < col.length; i++) {
      let colname = col[i]
      attrStack[colname] = i
    }
    var dataStack = src.map((srcitem) => (Object.fromEntries(col.map((item,index)=>([item,srcitem[index]])))))
    setData(dataStack)
    //filter columns
    let filterName = ['species', 'tissues', 'technologies']
    let filterItem = filterName.map((fc) => {
      let Ori_col = Array.from(new Set(dataStack.map((val) => val[fc])))
      let Uni_col = []
      for (let col of Ori_col) {
        if (col !== null && typeof col == 'string') {
          Uni_col = [...Uni_col, ...col?.split(';')]
        }
      }
      Uni_col = Array.from(new Set(Uni_col)).sort()
      Uni_col = Uni_col.map((item) => ({ text: item, value: item }))
      return Uni_col
    })
    let filterCol = {}
    for (let i = 0; i < filterName.length; i++) {
      filterCol[filterName[i]] = filterItem[i]
    }

    //generate columns and filters
    let columns = [
      Table.SELECTION_COLUMN,
      {
        title: 'ID',
        dataIndex: 'id',
        key: 'id',
        width: '4%',
        sorter: (a, b) => a.key - b.key,
        sortDirections: ['descend', 'ascend'],
        fixed: 'left',
      },
      {
        title: 'Title',
        dataIndex: 'title',
        ellipsis: true,
        key: 'title',
        width: '35%',
        ...getColumnSearchProps('title'),
      },
      {
        title: 'Species',
        dataIndex: 'species',
        key: 'species',
        ellipsis: true,
        width: '15%',
        filters: filterCol['species'],
        filterMode: 'tree',
        filterSearch: true,
        onFilter: (value, record) => record.species.split(';').includes(value),
      },
      {
        title: 'Tissues',
        dataIndex: 'tissues',
        key: 'tissues',
        ellipsis: true,
        width: '15%',
        filters: filterCol['tissues'],
        filterMode: 'tree',
        filterSearch: true,
        onFilter: (value, record) => {
          return record.tissues !== null
            ? record.tissues.split(';').includes(value)
            : false
        },
      },
      {
        title: 'Technology',
        dataIndex: 'technologies',
        key: 'technologies',
        width: '15%',
        ellipsis: true,
        filters: filterCol['technologies'],
        filterMode: 'tree',
        filterSearch: true,
        onFilter: (value, record) =>
          record.technologies.split(';').includes(value),
      },
      {
        title: 'PMID',
        dataIndex: 'pmid',
        ellipsis: true,
        key: 'pmid',
        ...getColumnSearchProps('pmid'),
        width: '10%',
      },
      {
        title: 'Action',
        key: 'operation',
        fixed: 'right',
        width: '10%',
        render: (index, record) => (
          <>
          {/* <div>{JSON.stringify(record.id)}</div> */}
            <Space>
              <DatasetDescription
                descCol={col}
                descInfo={Object.values(record)}
                placement={'top'}
              />
              <IconTip
                icon={EyeOutlined}
                text="Visualization"
                key="list-vertical-visualization-o"
                attr="Visualize this Dataset"
                placement="top"
                onClick={() => {
                  let scid = record['has_paired']
                  let scitem = GetscItem(scid)
                  let state = {
                    st: record
                  }
                  if (typeof scitem !== 'undefined') {
                    let scvalues = Object.fromEntries(
                      dataCol.map((k, i) => [k, scitem[i]])
                    )
                    state['sc'] = scvalues
                  }
                  // console.log(state)
                  navigate('/browse', { state: state })
                }}
              />
              <DatasetPaired
                descCol={col}
                descInfo={Object.values(record)}
                placement={'top'}
                // scInfo={item}
                scInfo={
                  src.find(s => s[1] === record['has_paired'])
                  // src.find(s => s[1] === item[26])
                  // Array.isArray(src) && item.length > 26
                  //   ? src.find(s => Array.isArray(s) && s[1] === item[26]) || item
                  //   : null
                }
              />
              <IconTip
                icon={CloudDownloadOutlined}
                attr={'Download datasets'}
                placement="top"
                onClick={() => {
                  // console.log(src)
                  // console.log(item[26])
                  let x = src.find(s => s[1] === record['has_paired'])
                  // console.log(x)
                  let url = record['accessions']  // 获取 item[23] 作为 URL
                  // console.log(url);
                  if (url) {
                    /*** 直接跳转数据集的标签页 ***/
                    window.open(url, '_blank')  // 在新标签页中打开链接
                  } else {
                    /*** 从我们的resource中直接下载 ***/
                    console.error('URL is invalid')
                    axios
                      .get('/api/query', {
                        responseType: 'blob',
                      })
                      .then((response) => {
                        // console.log(response)
                        let blob = response.data
                        let url = window.URL.createObjectURL(blob)
                        let a = document.createElement('a')
                        a.href = url
                        a.download = `${record['dataset_id']}.h5ad`
                        document.body.appendChild(a)
                        a.click()

                        // 清理
                        window.URL.revokeObjectURL(url)
                        document.body.removeChild(a)
                      })
                      .catch(error => {
                        console.error('Error fetching blob:', error)
                      })
                  }
                }}
              />
            </Space>
          </>
        ),
      },
    ]
    setDyColumns(columns)
  }, [src, col])

  return (
    <>
      <Table
        columns={dyColumns}
        dataSource={data}
        size="small"
        rowSelection={{}}
        scroll={{
          x: 1000,
          y: 1000,
        }}
        pagination={{
          showSizeChanger: true,
          showQuickJumper: true,
          showTotal: (total) => `Total ${total} items`,
          position: 'bottomRight',
        }}
      />
    </>
  )
}
export default DatasetTable
