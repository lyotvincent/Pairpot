import {
  InboxOutlined,
  DeleteOutlined,
  CheckCircleOutlined,
} from '@ant-design/icons'
import { message, Upload } from 'antd'
import React, { useState, useEffect } from 'react'
import styles from './upload.less'
import { Stack, Card, Toast, Row, Col } from 'react-bootstrap'
import { LoadingButton } from './LoadingButton'
const { Dragger } = Upload

export const FileUpload = () => {
  const [fileList, setFileList] = useState([])
  const [placeholder, setPlaceholder] = useState(
    'Support .h5spt, .h5ad or .zip file with 10X Visium format.'
  )

  const handleChange = (info) => {
    var newFileList = [...info.fileList] // load the previous fileList

    console.log(info)

    newFileList = newFileList.map((file) => {
      //Read from response and show file link
      if (file.response) {
        // Component will show file.url as link
        file.url = file.response.url
      }
      return file
    })
    setFileList(newFileList)
  }

  const handleDrop = (e) => {
    console.log('Dropped files', e.dataTransfer.files)
  }

  const props = {
    name: 'file',
    action: 'https://www.mocky.io/v2/5cc8019d300000980a055e76',
    onChange: handleChange,
    onDrop: handleDrop,
    multiple: false,
    showUploadList: false,
  }

  useEffect(() => {
    console.log('placeholder: %s', placeholder)
    var newPlaceholder = fileList.map((file) => {
      return (
        <Toast.Body>
          <Stack direction="horizontal" gap={3}>
            {' '}
            <CheckCircleOutlined />
            {file.name} <DeleteOutlined />
          </Stack>
        </Toast.Body>
      )
    })
    setPlaceholder(newPlaceholder)
  }, fileList)

  return (
    <div style={{ margin: '2rem' }}>
      <Row>
        <Col>
          <Dragger {...props} fileList={fileList}>
            <p className="ant-upload-drag-icon">
              <InboxOutlined color="red" />
            </p>
            <p>Click or drag file to this area to upload</p>
            <p>
              Support for a single or bulk upload. Strictly prohibit from
              uploading company data or other band files
            </p>
          </Dragger>
        </Col>
        <Col>
          <Toast bg="selected">
            <Toast.Header closeButton={true}>
              <strong className="me-auto">dataset</strong>
              <small>1 min ago</small>
            </Toast.Header>
            {placeholder}
          </Toast>
        </Col>
      </Row>

      <br />
      <LoadingButton value="Browse..." variant="success" />
    </div>
  )
}
