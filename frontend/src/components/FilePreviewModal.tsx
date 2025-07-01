import React, { useEffect, useState, useRef } from 'react';
import Modal from 'react-modal';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

Modal.setAppElement('#root'); // Adjust if your app root is different

interface FilePreviewModalProps {
  fileId: string | null;
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
      fetch(`/api/files/${fileId}/inspect`, { credentials: 'include' })
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
        return <pre style={{ whiteSpace: 'pre-wrap', maxHeight: 650, overflowY: 'auto' }}>{fileData.content}</pre>;
      case 'markdown':
        return (
          <div
            style={{
              maxHeight: 650,
              overflowY: 'auto',
              lineHeight: 1.5,
              fontSize: 16,
              padding: 16,
              background: '#20232a',
              borderRadius: 8,
            }}
            className="markdown-body"
          >
            <ReactMarkdown remarkPlugins={[remarkGfm]}>{fileData.content}</ReactMarkdown>
          </div>
        );
      case 'csv':
        return (
          <div ref={tableScrollRef} style={{ maxHeight: 650, overflowY: 'auto', paddingLeft: 32, paddingRight: 32 }}>
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
        return <pre style={{ maxHeight: 650, overflowY: 'auto' }}>{JSON.stringify(fileData.content, null, 2)}</pre>;
      case 'xlsx':
        return (
          <div ref={tableScrollRef} style={{ maxHeight: 650, overflowY: 'auto', paddingLeft: 32 }}>
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
        return <img src={fileData.content} alt={fileData.filename} style={{ maxWidth: '100%', maxHeight: 220 }} />;
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
        content: { maxWidth: 1400, margin: 'auto', minHeight: 400, maxHeight: 800, background: '#181f1b', color: '#fff', borderRadius: 12, overflow: 'auto', top: '5%' },
        overlay: { backgroundColor: 'rgba(0,0,0,0.7)' }
      }}
    >
      <h2 className="text-lg font-bold mb-4">File Preview</h2>
      {fileData?.filename && (
        <div className="mb-2 text-sm text-gray-300 font-mono truncate flex items-center gap-4" title={fileData.filename}>
          <span>{fileData.filename}</span>
          <a
            href={`/api/files/${fileId}/download`}
            className="px-3 py-1 border border-gray-400 text-gray-300 bg-transparent hover:bg-green-400/10 hover:border-green-300 hover:text-green-200 rounded-full text-xs font-semibold flex items-center gap-1 transition focus:outline-none focus:ring-2 focus:ring-green-400/50"
            title={`Download ${fileData.filename}`}
            download
            style={{ textDecoration: 'none' }}
          >
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 20 20" fill="currentColor" className="w-4 h-4 mr-1"><path d="M9 3.75a1 1 0 0 1 2 0v6.19l1.72-1.72a1 1 0 1 1 1.42 1.42l-3.43 3.43a1 1 0 0 1-1.42 0l-3.43-3.43a1 1 0 1 1 1.42-1.42L9 9.94V3.75ZM4.25 15a.75.75 0 0 1 .75-.75h10a.75.75 0 0 1 0 1.5H5a.75.75 0 0 1-.75-.75Z" /></svg>
            Download
          </a>
        </div>
      )}
      {loading && <div>Loading...</div>}
      {error && <div style={{ color: 'red' }}>{error}</div>}
      {!loading && !error && renderPreview()}
    </Modal>
  );
};

export default FilePreviewModal; 