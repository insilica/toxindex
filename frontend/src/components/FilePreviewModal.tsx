import React, { useEffect, useState, useRef } from 'react';
import Modal from 'react-modal';

Modal.setAppElement('#root'); // Adjust if your app root is different

interface FilePreviewModalProps {
  fileId: number | null;
  isOpen: boolean;
  onRequestClose: () => void;
}

const FilePreviewModal: React.FC<FilePreviewModalProps> = ({ fileId, isOpen, onRequestClose }) => {
  const [fileData, setFileData] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const tableScrollRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (isOpen && fileId) {
      setLoading(true);
      setError('');
      fetch(`/api/file/${fileId}/inspect`, { credentials: 'include' })
        .then(res => res.json())
        .then(data => {
          setFileData(data);
          setLoading(false);
        })
        .catch(() => {
          setError('Failed to load file preview');
          setLoading(false);
        });
    } else if (!isOpen) {
      setFileData(null);
      setError('');
    }
  }, [isOpen, fileId]);

  const renderPreview = () => {
    if (!fileData) return null;
    switch (fileData.type) {
      case 'text':
        return <pre style={{ whiteSpace: 'pre-wrap', maxHeight: 400, overflowY: 'auto' }}>{fileData.content}</pre>;
      case 'markdown':
        return <div style={{ maxHeight: 400, overflowY: 'auto' }} dangerouslySetInnerHTML={{ __html: fileData.content }} />;
      case 'csv':
        return (
          <div ref={tableScrollRef} style={{ maxHeight: 400, overflowY: 'auto', paddingLeft: 32, paddingRight: 32 }}>
            <table className="min-w-full text-xs text-gray-200">
              <thead>
                <tr>
                  {fileData.content[0]?.map((cell: string, j: number) => (
                    <th
                      key={j}
                      className="border px-2 py-1 bg-gray-800 sticky top-0 z-10 truncate"
                      style={{
                        background: '#23272b',
                        position: 'sticky',
                        top: 0,
                        maxWidth: 120,
                        whiteSpace: 'nowrap',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        borderTop: '1px solid #fff',
                        boxShadow: '0 1px 0 0 #fff'
                      }}
                    >
                      {cell}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {fileData.content.slice(1).map((row: string[], i: number) => (
                  <tr key={i}>
                    {row.map((cell: string, j: number) => (
                      <td
                        key={j}
                        className="border px-2 py-1 truncate"
                        style={{ maxWidth: 120, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}
                      >
                        {cell}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        );
      case 'json':
        return <pre style={{ maxHeight: 400, overflowY: 'auto' }}>{JSON.stringify(fileData.content, null, 2)}</pre>;
      case 'xlsx':
        return (
          <div ref={tableScrollRef} style={{ maxHeight: 400, overflowY: 'auto', paddingLeft: 32 }}>
            <table className="min-w-full text-xs text-gray-200">
              <thead>
                <tr>
                  {fileData.columns.map((col: string, i: number) => (
                    <th
                      key={i}
                      className="border px-2 py-1 bg-gray-800 sticky top-0 z-10 truncate"
                      style={{
                        background: '#23272b',
                        position: 'sticky',
                        top: 0,
                        maxWidth: 120,
                        whiteSpace: 'nowrap',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                        boxShadow: '0 3px 0 0 #fff'
                      }}
                    >
                      {col}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {fileData.content.map((row: any, i: number) => (
                  <tr key={i}>
                    {fileData.columns.map((col: string, j: number) => (
                      <td
                        key={j}
                        className="border px-2 py-1 truncate"
                        style={{ maxWidth: 120, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}
                      >
                        {row[col]}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        );
      case 'image':
        return <img src={fileData.content} alt={fileData.filename} style={{ maxWidth: '100%', maxHeight: 400 }} />;
      case 'unsupported':
        return <div>Preview not supported for this file type.</div>;
      case 'error':
        return <div>Error: {fileData.error}</div>;
      default:
        return <div>Unknown file type.</div>;
    }
  };

  return (
    <Modal
      isOpen={isOpen}
      onRequestClose={onRequestClose}
      contentLabel="File Preview"
      style={{
        content: { maxWidth: 800, margin: 'auto', minHeight: 200, background: '#181f1b', color: '#fff', borderRadius: 12 },
        overlay: { backgroundColor: 'rgba(0,0,0,0.7)' }
      }}
    >
      <h2 className="text-lg font-bold mb-4">File Preview</h2>
      {loading && <div>Loading...</div>}
      {error && <div style={{ color: 'red' }}>{error}</div>}
      {!loading && !error && renderPreview()}
    </Modal>
  );
};

export default FilePreviewModal; 