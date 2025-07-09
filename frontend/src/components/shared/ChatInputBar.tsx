import React, { useState } from 'react';
import { FaPlus } from 'react-icons/fa';
import EnvironmentSelector from './EnvironmentSelector';
import UploadCsvModal from './UploadCsvModal';
import { useEnvironment } from "../../context/EnvironmentContext";
import WorkflowSelector from './WorkflowSelector';
import { FolderSearch } from 'lucide-react';
import FilePickerModal from './FilePickerModal';

interface ChatInputBarProps {
  value: string;
  onChange: (v: string) => void;
  onSubmit: (e: React.FormEvent<HTMLFormElement>) => void;
  uploading?: boolean;
  error?: string | null;
  placeholder?: string;
  disabled?: boolean;
  onFilePick?: (fileId: string, fileName?: string) => void;
}

const ChatInputBar: React.FC<ChatInputBarProps> = ({
  value,
  onChange,
  onSubmit,
  uploading = false,
  error,
  placeholder = 'Ask me your toxicology question [ Is green tea nephrotoxic? ]',
  onFilePick,
}) => {
  const [showUploadModal, setShowUploadModal] = useState(false);
  const [showFilePicker, setShowFilePicker] = useState(false);
  const { environments, selectedEnv } = useEnvironment();

  const handleUploadClick = () => {
    setShowUploadModal(true);
  };

  return (
    <>
      <div className="w-full z-9 bg-opacity-10 flex justify-center" style={{ padding: 0 }}>
        <form className="flex flex-col items-center w-full" style={{ maxWidth: '800px' }} onSubmit={onSubmit}>
          <div className="relative w-full" style={{ maxWidth: '800px', width: '100%' }}>
            <textarea
              rows={3}
              placeholder={placeholder}
              value={value}
              onChange={e => onChange(e.target.value)}
              onKeyDown={e => {
                if (e.key === 'Enter' && !e.shiftKey) {
                  e.preventDefault();
                  onSubmit(e as any);
                }
              }}
              className="w-full pt-4 pb-4 pr-28 pl-8 text-lg rounded-2xl border border-gray-700 bg-gray-900 bg-opacity-70 text-white resize-none min-h-[80px] shadow-2xl focus:outline-none focus:ring-2 focus:ring-green-400 text-left placeholder:text-left"
              style={{ minHeight: 80, fontFamily: 'inherit', width: '100%', boxShadow: '0 8px 32px 0 rgba(34,197,94,0.10)' }}
              disabled={uploading}
            />
            <div className="absolute right-4 bottom-14 flex items-center space-x-2 z-10">
              <button
                className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                title="Pick a file from environment"
                onClick={e => { e.preventDefault(); setShowFilePicker(true); }}
                disabled={uploading}
                type="button"
              >
                <FolderSearch className="w-5 h-5 text-black" />
              </button>
              <button
                className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                title="Upload CSV file"
                onClick={e => { e.preventDefault(); handleUploadClick(); }}
                disabled={uploading}
                type="button"
              >
                <FaPlus className="w-5 h-5 text-black" />
              </button>
              <button
                type="submit"
                className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                title="Submit"
                disabled={uploading}
              >
                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={2.2} stroke="black" className="w-5 h-5">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M5 10l7-7m0 0l7 7m-7-7v18" />
                </svg>
              </button>
            </div>
            <div className="absolute left-4 z-10 flex items-end gap-2" style={{ position: 'relative', width: '420px', overflow: 'visible', bottom: '3rem' }}>
              <EnvironmentSelector/>
              <WorkflowSelector />
            </div>
          </div>
          {error && <div className="text-red-500 mt-2">{error}</div>}
        </form>
      </div>
      <UploadCsvModal
        open={showUploadModal}
        onClose={() => setShowUploadModal(false)}
        environments={environments}
        defaultEnvId={selectedEnv || undefined}
        onUploadSuccess={() => setShowUploadModal(false)}
      />
      <FilePickerModal
        open={showFilePicker}
        onClose={() => setShowFilePicker(false)}
        environmentId={selectedEnv}
        environments={environments}
        onPick={(fileId: string, fileName?: string) => {
          setShowFilePicker(false);
          onFilePick?.(fileId, fileName);
        }}
      />
    </>
  );
};

export default ChatInputBar; 