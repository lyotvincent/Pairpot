import React, { useEffect, useState } from 'react'
import { Button, Stack } from 'react-bootstrap'
import { LoadingOutlined } from '@ant-design/icons'

function simulateNetworkRequest() {
  return new Promise((resolve) => setTimeout(resolve, 2000))
}

export const LoadingButton = (props) => {
  const [loading, setLoading] = useState(false)

  useEffect(() => {
    if (loading) {
      simulateNetworkRequest().then(() => {
        setLoading(false)
      })
    }
  }, [loading])

  const handleClick = () => setLoading(true)

  var loadingElement = (
    <Stack direction="horizontal" gap={2}>
      <LoadingOutlined /> loading...
    </Stack>
  )
  return (
    <Button
      variant={props.variant}
      disabled={loading}
      onClick={loading ? null : handleClick}>
      {loading ? loadingElement : props.value}
    </Button>
  )
}
